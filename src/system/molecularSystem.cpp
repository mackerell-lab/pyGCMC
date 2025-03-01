#include "system/molecularSystem.hpp"
#include "system/system.hpp"
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <map>
#include <set>
#include <queue>

namespace pygcmc {
namespace system {

using System = pygcmc::system::System;
using LogLevel = pygcmc::system::LogLevel;

std::shared_ptr<model::Molecular> MolecularSystem::combine(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology) {
    
    if (!structure || !topology) {
        throw std::invalid_argument("Structure and Topology pointers cannot be null");
    }

    // Create new Molecular object
    molecular_ = std::make_shared<model::Molecular>();

    // Copy data from Structure
    molecular_->atoms = structure->get_atoms();
    molecular_->residues = structure->get_residues();
    molecular_->terminals = structure->get_terminals();
    molecular_->helices = structure->get_helices();
    molecular_->sheets = structure->get_sheets();
    molecular_->ssbonds = structure->get_ssbonds();
    molecular_->boxDimensions = structure->get_box_dimensions();

    // Copy data from Topology
    const auto num_atoms = static_cast<size_t>(topology->get_num_atoms());
    const auto num_residues = static_cast<size_t>(topology->get_num_residues());
    const auto num_segments = static_cast<size_t>(topology->get_num_segments());

    // Verify total atom count
    if (molecular_->atoms.size() != num_atoms) {
        std::stringstream ss;
        ss << "Inconsistent total number of atoms: Structure has " 
           << molecular_->atoms.size() << " atoms, but Topology has " 
           << num_atoms << " atoms\n"
           << "This mismatch suggests that the structure file contains additional molecules "
           << "that are not present in the topology file.\n"
           << "Structure residues (in PDB order):";
        
        // List residues from structure file in PDB order
        for (const auto& res : molecular_->residues) {
            ss << "\n  " << res->get_resname() << " " << res->get_ires();
        }
        
        ss << "\n\nTopology residues:";
        // List residues from topology file in order
        for (size_t i = 0; i < num_residues; ++i) {
            const auto& res = topology->get_residue(static_cast<int>(i));
            ss << "\n  " << res.name << " " << res.number;
        }
        
        // Add residue statistics
        ss << "\n\nResidue count summary:";
        ss << "\nStructure:";
        std::map<std::string, int> struct_res_count;
        for (const auto& res : molecular_->residues) {
            struct_res_count[res->get_resname()]++;
        }
        for (const auto& [resname, count] : struct_res_count) {
            ss << "\n  " << resname << ": " << count;
        }
        
        ss << "\nTopology:";
        std::map<std::string, int> top_res_count;
        for (size_t i = 0; i < num_residues; ++i) {
            const auto& res = topology->get_residue(static_cast<int>(i));
            top_res_count[res.name]++;
        }
        for (const auto& [resname, count] : top_res_count) {
            ss << "\n  " << resname << ": " << count;
        }
        
        throw std::runtime_error(ss.str());
    }

    // Verify total residue count
    if (molecular_->residues.size() != num_residues) {
        std::stringstream ss;
        ss << "Inconsistent total number of residues: Structure has " 
           << molecular_->residues.size() << " residues, but Topology has " 
           << num_residues << " residues\n"
           << "This mismatch suggests that the structure file contains additional residues "
           << "that are not present in the topology file.\n"
           << "Structure residues (in PDB order):";
        
        // List residues from structure file in PDB order
        for (const auto& res : molecular_->residues) {
            ss << "\n  " << res->get_resname() << " " << res->get_ires();
        }
        
        ss << "\n\nTopology residues:";
        // List residues from topology file in order
        for (size_t i = 0; i < num_residues; ++i) {
            const auto& res = topology->get_residue(static_cast<int>(i));
            ss << "\n  " << res.name << " " << res.number;
        }
        
        // Add residue statistics
        ss << "\n\nResidue count summary:";
        ss << "\nStructure:";
        std::map<std::string, int> struct_res_count;
        for (const auto& res : molecular_->residues) {
            struct_res_count[res->get_resname()]++;
        }
        for (const auto& [resname, count] : struct_res_count) {
            ss << "\n  " << resname << ": " << count;
        }
        
        ss << "\nTopology:";
        std::map<std::string, int> top_res_count;
        for (size_t i = 0; i < num_residues; ++i) {
            const auto& res = topology->get_residue(static_cast<int>(i));
            top_res_count[res.name]++;
        }
        for (const auto& [resname, count] : top_res_count) {
            ss << "\n  " << resname << ": " << count;
        }
        
        throw std::runtime_error(ss.str());
    }

    // Verify atom type first letter match
    for (size_t i = 0; i < num_atoms; ++i) {
        const auto& pdb_atom = molecular_->atoms[i];
        const auto& top_atom = topology->get_atom(static_cast<int>(i));
        
        // Extract element from PDB atom name
        std::string pdb_element = pdb_atom->get_element();
        if (pdb_element.empty()) {
            pdb_element = pdb_atom->get_type();
        }
        
        // Extract element from topology atom type
        std::string top_element = top_atom.type;
        
        // Compare first letter (converted to uppercase)
        char pdb_first = std::toupper(pdb_element[0]);
        char top_first = std::toupper(top_element[0]);
        
        if (pdb_first != top_first) {
            std::stringstream ss;
            ss << "Mismatched atom elements at index " << i << ": Structure has " 
               << pdb_first << " (from " << pdb_atom->get_type() 
               << "), but Topology has " << top_first 
               << " (from " << top_atom.type << ")";
            throw std::runtime_error(ss.str());
        }
    }

    // Standard amino acid list
    static const std::set<std::string> standard_amino_acids = {
        "ALA", "ARG", "ASN", "ASP", "CYS", 
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
        "HSE", "HSP", "HSC",  // Histidine different protonation states
        "CYX",  // Disulfide bond form cysteine
        "HID", "HIE", "HIP"   // Histidine in CHARMM force field
    };

    // Use bond information to build residue connection relationships
    std::map<int, std::set<int>> residue_connections;  // residue_id -> connected residue_ids
    const auto& bonds = topology->get_bonds();
    for (const auto& bond : bonds) {
        const auto& atom1 = topology->get_atom(bond.atom1);
        const auto& atom2 = topology->get_atom(bond.atom2);
        if (atom1.residue_id != atom2.residue_id) {
            residue_connections[atom1.residue_id].insert(atom2.residue_id);
            residue_connections[atom2.residue_id].insert(atom1.residue_id);
        }
    }

    // Use BFS to find all connected residue groups (chains)
    std::set<int> visited_residues;
    std::vector<std::set<int>> chains;  // Each element is a list of residue IDs in a chain
    
    for (size_t i = 0; i < num_residues; ++i) {
        int res_id = static_cast<int>(i);
        if (visited_residues.find(res_id) != visited_residues.end()) {
            continue;
        }

        // Find all residues connected to the current residue
        std::set<int> current_chain;
        std::queue<int> to_visit;
        to_visit.push(res_id);
        
        while (!to_visit.empty()) {
            int current = to_visit.front();
            to_visit.pop();
            
            if (visited_residues.find(current) != visited_residues.end()) {
                continue;
            }
            
            visited_residues.insert(current);
            current_chain.insert(current);
            
            // Add connected residues to the queue
            if (residue_connections.find(current) != residue_connections.end()) {
                for (int connected : residue_connections[current]) {
                    if (visited_residues.find(connected) == visited_residues.end()) {
                        to_visit.push(connected);
                    }
                }
            }
        }
        
        if (!current_chain.empty()) {
            chains.push_back(current_chain);
        }
    }

    // Calculate amino acid count on each chain
    std::map<std::set<int>, int> chain_aa_count;  // chain -> amino acid count
    for (const auto& chain : chains) {
        int aa_count = 0;
        for (int res_id : chain) {
            const auto& res = topology->get_residue(res_id);
            if (standard_amino_acids.find(res.name) != standard_amino_acids.end()) {
                aa_count++;
            }
        }
        chain_aa_count[chain] = aa_count;
    }

    // Establish mapping for non-protein residues
    std::map<std::string, std::vector<size_t>> mol_type_indices;  // residue name -> list of indices in topology
    for (size_t i = 0; i < num_residues; ++i) {
        const auto& top_res = topology->get_residue(static_cast<int>(i));
        mol_type_indices[top_res.name].push_back(i);
    }

    // Maintain current index for each type of non-protein molecule
    std::map<std::string, size_t> current_mol_index;

    // Verify atom count and type for each residue
    for (const auto& res : molecular_->residues) {
        const auto& res_atoms = res->get_atoms();
        bool found_matching_res = false;
        
        // Find current residue's chain
        int current_res_id = -1;
        std::set<int>* current_chain = nullptr;
        
        for (size_t i = 0; i < num_residues; ++i) {
            const auto& top_res = topology->get_residue(static_cast<int>(i));
            if (top_res.name == res->get_resname()) {
                current_res_id = static_cast<int>(i);
                break;
            }
        }
        
        if (current_res_id >= 0) {
            for (auto& chain : chains) {
                if (chain.find(current_res_id) != chain.end()) {
                    current_chain = &chain;
                    break;
                }
            }
        }
        
        // Check if it's a protein residue
        bool is_amino_acid = standard_amino_acids.find(res->get_resname()) != standard_amino_acids.end();
        bool is_protein = false;
        if (current_chain && chain_aa_count[*current_chain] > 3) {
            is_protein = is_amino_acid;
        }
        
        if (is_protein) {
            // For protein residues, match by name and number exactly
            for (size_t i = 0; i < num_residues; ++i) {
                const auto& top_res = topology->get_residue(static_cast<int>(i));
                if (top_res.name == res->get_resname() && top_res.number == res->get_ires()) {
                    found_matching_res = true;
                    // Verify atom count
                    if (res_atoms.size() != top_res.atoms.size()) {
                        std::stringstream ss;
                        ss << "Inconsistent number of atoms in residue " << res->get_resname() 
                           << " " << res->get_ires() << ": Structure has " 
                           << res_atoms.size() << " atoms, but Topology has " 
                           << top_res.atoms.size() << " atoms";
                        throw std::runtime_error(ss.str());
                    }
                    // Verify atom type
                    for (size_t j = 0; j < res_atoms.size(); ++j) {
                        const auto& pdb_atom = res_atoms[j];
                        const auto& top_atom = topology->get_atom(top_res.atoms[j]);
                        
                        // Extract element from PDB atom name
                        std::string pdb_element = pdb_atom->get_element();
                        if (pdb_element.empty()) {
                            pdb_element = pdb_atom->get_type();
                        }
                        
                        // Extract element from topology atom type
                        std::string top_element = top_atom.type;
                        
                        // Compare first letter (converted to uppercase)
                        char pdb_first = std::toupper(pdb_element[0]);
                        char top_first = std::toupper(top_element[0]);
                        
                        if (pdb_first != top_first) {
                            std::stringstream ss;
                            ss << "Mismatched atom elements in residue " << res->get_resname() 
                               << " " << res->get_ires() << ": Structure has " 
                               << pdb_first << " (from " << pdb_atom->get_type() 
                               << "), but Topology has " << top_first 
                               << " (from " << top_atom.type << ")";
                            throw std::runtime_error(ss.str());
                        }
                    }
                    break;
                }
            }
        } else {
            // For non-protein residues, match by molecule type in sequence
            const std::string& resname = res->get_resname();
            auto it = mol_type_indices.find(resname);
            if (it != mol_type_indices.end()) {
                System::log(LogLevel::DEBUG, "Found ", it->second.size(), " instances of ", 
                          resname, " in topology, current index: ", current_mol_index[resname]);
            } else {
                System::log(LogLevel::DEBUG, "No instances of ", resname, " found in topology");
            }
            
            if (it != mol_type_indices.end() && current_mol_index[resname] < it->second.size()) {
                size_t top_res_idx = it->second[current_mol_index[resname]];
                const auto& top_res = topology->get_residue(static_cast<int>(top_res_idx));
                found_matching_res = true;
                current_mol_index[resname]++;
                
                // Verify atom count
                if (res_atoms.size() != top_res.atoms.size()) {
                    std::stringstream ss;
                    ss << "Inconsistent number of atoms in residue " << res->get_resname() 
                       << " " << res->get_ires() << ": Structure has " 
                       << res_atoms.size() << " atoms, but Topology has " 
                       << top_res.atoms.size() << " atoms";
                    throw std::runtime_error(ss.str());
                }
                // Check atom types
                for (size_t j = 0; j < res_atoms.size(); ++j) {
                    const auto& pdb_atom = res_atoms[j];
                    const auto& top_atom = topology->get_atom(top_res.atoms[j]);
                    
                    // Extract element from PDB atom name
                    std::string pdb_element = pdb_atom->get_element();
                    if (pdb_element.empty()) {
                        pdb_element = pdb_atom->get_type();
                    }
                    
                    // Extract element from topology atom type
                    std::string top_element = top_atom.type;
                    
                    // Compare first letter (converted to uppercase)
                    char pdb_first = std::toupper(pdb_element[0]);
                    char top_first = std::toupper(top_element[0]);
                    
                    if (pdb_first != top_first) {
                        std::stringstream ss;
                        ss << "Mismatched atom elements in residue " << res->get_resname() 
                           << " " << res->get_ires() << ": Structure has " 
                           << pdb_first << " (from " << pdb_atom->get_type() 
                           << "), but Topology has " << top_first 
                           << " (from " << top_atom.type << ")";
                        throw std::runtime_error(ss.str());
                    }
                }
            }
        }
        
        if (!found_matching_res) {
            std::stringstream ss;
            ss << "Could not find matching residue in topology for " 
               << res->get_resname() << " " << res->get_ires();
            throw std::runtime_error(ss.str());
        }
    }

    molecular_->topology_atoms.reserve(num_atoms);
    molecular_->topology_residues.reserve(num_residues);
    molecular_->segments.reserve(num_segments);

    for (size_t i = 0; i < num_atoms; ++i) {
        molecular_->topology_atoms.push_back(topology->get_atom(static_cast<int>(i)));
    }
    for (size_t i = 0; i < num_residues; ++i) {
        molecular_->topology_residues.push_back(topology->get_residue(static_cast<int>(i)));
    }
    for (size_t i = 0; i < num_segments; ++i) {
        molecular_->segments.push_back(topology->get_segment(static_cast<int>(i)));
    }

    // Copy bonding information
    molecular_->bonds = topology->get_bonds();
    molecular_->angles = topology->get_angles();
    molecular_->dihedrals = topology->get_dihedrals();
    molecular_->donors = topology->get_donors();
    molecular_->acceptors = topology->get_acceptors();
    molecular_->exclusions = topology->get_exclusions();
    molecular_->groups = topology->get_groups();
    molecular_->cmaps = topology->get_cmaps();
    // Add standardized CMAPs
    for (const auto& cmap : topology->get_cmaps()) {
        molecular_->add_standard_cmap(cmap);
    }
    molecular_->titles = topology->get_titles();

    // Copy lookup mappings
    for (size_t i = 0; i < num_segments; ++i) {
        const auto& segment = topology->get_segment(static_cast<int>(i));
        molecular_->segment_map[segment.name] = segment.id;
    }

    for (size_t i = 0; i < num_residues; ++i) {
        const auto& residue = topology->get_residue(static_cast<int>(i));
        molecular_->residue_map[std::make_pair(residue.name, residue.number)] = residue.id;
    }

    for (size_t i = 0; i < num_atoms; ++i) {
        const auto& atom = topology->get_atom(static_cast<int>(i));
        const auto& residue = topology->get_residue(atom.residue_id);
        molecular_->atom_map[std::make_tuple(residue.name, residue.number, atom.name)] = atom.id;
    }

    return molecular_;
}

std::shared_ptr<model::Molecular> MolecularSystem::combine_multiple(
    const std::shared_ptr<model::Structure>& structure,
    const std::vector<std::shared_ptr<model::Topology>>& topologies) {
    
    if (!structure || topologies.empty()) {
        throw std::invalid_argument("Structure and Topologies cannot be null/empty");
    }

    // Create new Molecular object
    molecular_ = std::make_shared<model::Molecular>();

    // Copy data from Structure
    molecular_->atoms = structure->get_atoms();
    molecular_->residues = structure->get_residues();
    molecular_->terminals = structure->get_terminals();
    molecular_->helices = structure->get_helices();
    molecular_->sheets = structure->get_sheets();
    molecular_->ssbonds = structure->get_ssbonds();
    molecular_->boxDimensions = structure->get_box_dimensions();

    // Try to match each topology
    std::vector<bool> residue_matched(molecular_->residues.size(), false);
    std::vector<std::shared_ptr<model::Topology>> matched_topologies;
    
    // First try to match the longest residue sequence
    for (size_t start_idx = 0; start_idx < molecular_->residues.size(); ++start_idx) {
        if (residue_matched[start_idx]) continue;
        
        for (const auto& topology : topologies) {
            size_t matched_count = 0;
            if (match_residue_sequence(molecular_->residues, topology, start_idx, matched_count)) {
                // Mark matched residues
                for (size_t i = 0; i < matched_count; ++i) {
                    residue_matched[start_idx + i] = true;
                }
                matched_topologies.push_back(topology);
                break;
            }
        }
    }
    
    // Check if all residues have been matched
    for (size_t i = 0; i < residue_matched.size(); ++i) {
        if (!residue_matched[i]) {
            std::stringstream ss;
            ss << "Could not find matching topology for residue " 
               << molecular_->residues[i]->get_resname()
               << " " << molecular_->residues[i]->get_ires();
            throw std::runtime_error(ss.str());
        }
    }
    
    // Merge all matched topologies
    merge_topologies(molecular_, matched_topologies);
    
    return molecular_;
}

bool MolecularSystem::match_residue_sequence(
    const std::vector<std::shared_ptr<model::Residue>>& pdb_residues,
    const std::shared_ptr<model::Topology>& topology,
    size_t start_idx,
    size_t& matched_count) {
    
    matched_count = 0;
    const size_t top_num_residues = topology->get_num_residues();
    
    // If the remaining residue count is insufficient, return false directly
    if (start_idx + top_num_residues > pdb_residues.size()) {
        return false;
    }
    
    // Check if the residue sequence matches
    for (size_t i = 0; i < top_num_residues; ++i) {
        const auto& pdb_res = pdb_residues[start_idx + i];
        const auto& top_res = topology->get_residue(static_cast<int>(i));
        
        if (pdb_res->get_resname() != top_res.name) {
            return false;
        }
        
        try {
            verify_atom_types(pdb_res, top_res, topology);
        } catch (const std::runtime_error&) {
            return false;
        }
        
        matched_count++;
    }
    
    return true;
}

void MolecularSystem::merge_topologies(
    std::shared_ptr<model::Molecular>& molecular,
    const std::vector<std::shared_ptr<model::Topology>>& topologies) {
    
    size_t total_atoms = 0;
    size_t total_residues = 0;
    
    // Calculate totals
    for (const auto& topology : topologies) {
        total_atoms += topology->get_num_atoms();
        total_residues += topology->get_num_residues();
    }
    
    // Verify if totals match
    if (molecular->atoms.size() != total_atoms) {
        std::stringstream ss;
        ss << "Total number of atoms mismatch: Structure has "
           << molecular->atoms.size() << " atoms, but combined topologies have "
           << total_atoms << " atoms";
        throw std::runtime_error(ss.str());
    }
    
    if (molecular->residues.size() != total_residues) {
        std::stringstream ss;
        ss << "Total number of residues mismatch: Structure has "
           << molecular->residues.size() << " residues, but combined topologies have "
           << total_residues << " residues";
        throw std::runtime_error(ss.str());
    }
    
    // Merge topology data
    molecular->topology_atoms.clear();
    molecular->topology_residues.clear();
    molecular->segments.clear();
    molecular->bonds.clear();
    molecular->angles.clear();
    molecular->dihedrals.clear();
    molecular->donors.clear();
    molecular->acceptors.clear();
    molecular->exclusions.clear();
    molecular->groups.clear();
    
    // Save existing CMAP
    std::vector<model::TopologyCmap> existing_cmaps = molecular->cmaps;
    molecular->cmaps.clear();
    
    size_t atom_offset = 0;
    size_t residue_offset = 0;
    
    for (const auto& topology : topologies) {
        // Copy atoms
        for (int i = 0; i < topology->get_num_atoms(); ++i) {
            auto atom = topology->get_atom(i);
            atom.id += atom_offset;
            atom.residue_id += residue_offset;
            molecular->topology_atoms.push_back(atom);
        }
        
        // Copy residues
        for (int i = 0; i < topology->get_num_residues(); ++i) {
            auto residue = topology->get_residue(i);
            residue.id += residue_offset;
            for (auto& atom_id : residue.atoms) {
                atom_id += atom_offset;
            }
            molecular->topology_residues.push_back(residue);
        }
        
        // Copy bonding information
        for (const auto& bond : topology->get_bonds()) {
            model::TopologyBond new_bond = bond;
            new_bond.atom1 += atom_offset;
            new_bond.atom2 += atom_offset;
            molecular->bonds.push_back(new_bond);
        }
        
        // Copy angle information
        for (const auto& angle : topology->get_angles()) {
            model::TopologyAngle new_angle = angle;
            new_angle.atom1 += atom_offset;
            new_angle.atom2 += atom_offset;
            new_angle.atom3 += atom_offset;
            molecular->angles.push_back(new_angle);
        }
        
        // Copy dihedral information
        for (const auto& dihedral : topology->get_dihedrals()) {
            model::TopologyDihedral new_dihedral = dihedral;
            new_dihedral.atom1 += atom_offset;
            new_dihedral.atom2 += atom_offset;
            new_dihedral.atom3 += atom_offset;
            new_dihedral.atom4 += atom_offset;
            molecular->dihedrals.push_back(new_dihedral);
        }
        
        // Copy CMAP information
        for (const auto& cmap : topology->get_cmaps()) {
            model::TopologyCmap new_cmap = cmap;
            
            // First update indices for all 8 atoms
            for (size_t i = 0; i < new_cmap.atoms.size(); ++i) {
                if (new_cmap.atoms[i] >= 0) {
                    const auto& atom = topology->get_atom(new_cmap.atoms[i]);
                    const auto& res = topology->get_residue(atom.residue_id);
                    if (res.name != "SOL") {
                        new_cmap.atoms[i] += atom_offset;
                    }
                }
            }
            
            molecular->cmaps.push_back(new_cmap);
            // Add standardized CMAP
            molecular->add_standard_cmap(new_cmap);
        }
        
        // Re-add previously saved CMAPs
        for (const auto& cmap : existing_cmaps) {
            molecular->cmaps.push_back(cmap);
            molecular->add_standard_cmap(cmap);
        }
        
        // Update offsets
        atom_offset += topology->get_num_atoms();
        residue_offset += topology->get_num_residues();
    }
}

void MolecularSystem::verify_atom_types(
    const std::shared_ptr<model::Residue>& pdb_res,
    const model::TopologyResidue& top_res,
    const std::shared_ptr<model::Topology>& topology) {
    
    const auto& pdb_atoms = pdb_res->get_atoms();
    if (pdb_atoms.size() != top_res.atoms.size()) {
        std::stringstream ss;
        ss << "Inconsistent number of atoms in residue " << pdb_res->get_resname()
           << " " << pdb_res->get_ires() << ": Structure has "
           << pdb_atoms.size() << " atoms, but Topology has "
           << top_res.atoms.size() << " atoms";
        throw std::runtime_error(ss.str());
    }
    
    for (size_t i = 0; i < pdb_atoms.size(); ++i) {
        const auto& pdb_atom = pdb_atoms[i];
        const auto& top_atom = topology->get_atom(static_cast<int>(top_res.atoms[i]));
        
        // Extract element from PDB atom name
        std::string pdb_element = pdb_atom->get_element();
        if (pdb_element.empty()) {
            pdb_element = pdb_atom->get_type();
        }
        
        // Extract element from topology atom type
        std::string top_element = top_atom.type;
        
        // Compare first letter (converted to uppercase)
        char pdb_first = std::toupper(pdb_element[0]);
        char top_first = std::toupper(top_element[0]);
        
        if (pdb_first != top_first) {
            std::stringstream ss;
            ss << "Mismatched atom elements in residue " << pdb_res->get_resname()
               << " " << pdb_res->get_ires() << ": Structure has "
               << pdb_first << " (from " << pdb_atom->get_type()
               << "), but Topology has " << top_first
               << " (from " << top_atom.type << ")";
            throw std::runtime_error(ss.str());
        }
    }
}

} // namespace system
} // namespace pygcmc
