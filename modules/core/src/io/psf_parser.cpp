// modules/core/src/io/psf_parser.cpp

#include "pygcmc/core/io/psf_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace pygcmc {
namespace core {
namespace io {

bool PSFParser::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::vector<std::string> atom_lines;
    bool in_atoms_section = false;
    std::string line;
    int natoms = 0;

    while (std::getline(file, line)) {
        // Skip empty lines and comments starting with *
        if (line.empty() || line[0] == '*') {
            continue;
        }

        // Check for NATOM section
        if (line.find("!NATOM") != std::string::npos) {
            // Extract number of atoms
            std::istringstream iss(line);
            iss >> natoms;
            if (natoms <= 0) {
                throw std::runtime_error("Invalid number of atoms in PSF file: " + std::to_string(natoms));
            }
            in_atoms_section = true;
            continue;
        }

        // End of atoms section
        if (in_atoms_section && line.find("!") != std::string::npos) {
            in_atoms_section = false;
            continue;
        }

        if (in_atoms_section && !line.empty()) {
            atom_lines.push_back(line);
        }
    }

    if (atom_lines.empty()) {
        throw std::runtime_error("No atoms found in PSF file");
    }

    return parse_atoms_section(atom_lines);
}

bool PSFParser::parse_atoms_section(const std::vector<std::string>& lines) {
    atoms_.clear();
    atom_index_.clear();
    id_index_.clear();

    for (const auto& line : lines) {
        std::istringstream iss(line);
        PSFAtom atom;

        // PSF format: ID SEGID RESID RESNAME ATOMNAME ATOMTYPE CHARGE MASS
        if (!(iss >> atom.id 
                  >> atom.segment 
                  >> atom.residue_number 
                  >> atom.residue 
                  >> atom.name 
                  >> atom.type 
                  >> atom.charge 
                  >> atom.mass)) 
        {
            std::cerr << "Failed to parse atom line: " << line << std::endl;
            continue;
        }

        std::cout << "Parsed atom: " << atom.residue << " " << atom.name 
                  << " type=" << atom.type 
                  << " charge=" << atom.charge 
                  << " mass=" << atom.mass << std::endl;

        // 存入容器
        size_t idx = atoms_.size();
        atoms_.push_back(atom);

        // 建立索引: residue_name -> (residue_number -> (atom_name -> idx))
        atom_index_[atom.residue][atom.residue_number][atom.name] = idx;
        
        // 建立 atom_id 索引
        id_index_[atom.id] = idx;
    }

    return !atoms_.empty();
}

bool PSFParser::get_atom_properties(const std::string& residue_name,
                                  const std::string& atom_name,
                                  double& charge,
                                  double& mass) const {
    // Find residue
    auto res_it = atom_index_.find(residue_name);
    if (res_it == atom_index_.end()) {
        std::cerr << "[get_atom_properties] Residue not found: " 
                  << residue_name << std::endl;
        return false;
    }

    // Get the first defined residue number for this residue type
    auto first_res_num_it = res_it->second.begin();
    if (first_res_num_it == res_it->second.end()) {
        std::cerr << "[get_atom_properties] No residue numbers found for: "
                  << residue_name << std::endl;
        return false;
    }

    // Find atom name
    auto atom_it = first_res_num_it->second.find(atom_name);
    if (atom_it == first_res_num_it->second.end()) {
        std::cerr << "[get_atom_properties] Atom not found: "
                  << residue_name << " " << atom_name << std::endl;
        return false;
    }

    const PSFAtom& atom = atoms_[atom_it->second];
    charge = atom.charge;
    mass = atom.mass;
    std::cout << "Found atom: " << residue_name << " " << atom_name 
              << " charge=" << charge << " mass=" << mass << std::endl;
    return true;
}

bool PSFParser::get_atom_properties(const std::string& residue_name,
                                  int residue_number,
                                  const std::string& atom_name,
                                  double& charge,
                                  double& mass) const {
    // Find residue
    auto res_it = atom_index_.find(residue_name);
    if (res_it == atom_index_.end()) {
        std::cerr << "[get_atom_properties] Residue not found: " 
                  << residue_name << std::endl;
        return false;
    }

    // Find residue number
    auto res_num_it = res_it->second.find(residue_number);
    if (res_num_it == res_it->second.end()) {
        std::cerr << "[get_atom_properties] Residue number not found: "
                  << residue_name << " " << residue_number << std::endl;
        return false;
    }

    // Find atom name
    auto atom_it = res_num_it->second.find(atom_name);
    if (atom_it == res_num_it->second.end()) {
        std::cerr << "[get_atom_properties] Atom not found: "
                  << residue_name << " " << residue_number << " " << atom_name << std::endl;
        return false;
    }

    const PSFAtom& atom = atoms_[atom_it->second];
    charge = atom.charge;
    mass = atom.mass;
    std::cout << "Found atom: " << residue_name << " " << residue_number << " " << atom_name 
              << " charge=" << charge << " mass=" << mass << std::endl;
    return true;
}

int PSFParser::update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const {
    int updated = 0;
    
    // First pass: match by exact residue number and atom name
    for (auto& pdb_atom : pdb_atoms) {
        for (const auto& psf_atom : atoms_) {
            if (psf_atom.residue == pdb_atom.residue && 
                psf_atom.residue_number == pdb_atom.sequence &&
                psf_atom.name == pdb_atom.name) {
                
                pdb_atom.topo_type = psf_atom.type;
                pdb_atom.topo_charge = psf_atom.charge;
                pdb_atom.topo_mass = psf_atom.mass;
                pdb_atom.chain = psf_atom.segment.empty() ? ' ' : psf_atom.segment[0];
                updated++;
                break;
            }
        }
    }
    
    // Second pass: match by residue name and atom name if sequence number didn't match
    for (auto& pdb_atom : pdb_atoms) {
        if (pdb_atom.topo_type.empty()) {  // Only try to match if not already matched
            for (const auto& psf_atom : atoms_) {
                if (psf_atom.residue == pdb_atom.residue && 
                    psf_atom.name == pdb_atom.name) {
                    
                    pdb_atom.topo_type = psf_atom.type;
                    pdb_atom.topo_charge = psf_atom.charge;
                    pdb_atom.topo_mass = psf_atom.mass;
                    pdb_atom.chain = psf_atom.segment.empty() ? ' ' : psf_atom.segment[0];
                    updated++;
                    break;
                }
            }
        }
    }

    return updated;
}

int PSFParser::update_pdb_atoms(std::vector<PDBAtom*>& pdb_atoms) const {
    int updated = 0;
    
    // First pass: match by exact residue number and atom name
    for (auto* pdb_atom : pdb_atoms) {
        if (!pdb_atom) continue;
        
        for (const auto& psf_atom : atoms_) {
            if (psf_atom.residue == pdb_atom->residue && 
                psf_atom.residue_number == pdb_atom->sequence &&
                psf_atom.name == pdb_atom->name) {
                
                pdb_atom->topo_type = psf_atom.type;
                pdb_atom->topo_charge = psf_atom.charge;
                pdb_atom->topo_mass = psf_atom.mass;
                pdb_atom->chain = psf_atom.segment.empty() ? ' ' : psf_atom.segment[0];
                updated++;
                break;
            }
        }
    }
    
    // Second pass: match by residue name and atom name if sequence number didn't match
    for (auto* pdb_atom : pdb_atoms) {
        if (!pdb_atom) continue;
        
        if (pdb_atom->topo_type.empty()) {  // Only try to match if not already matched
            for (const auto& psf_atom : atoms_) {
                if (psf_atom.residue == pdb_atom->residue && 
                    psf_atom.name == pdb_atom->name) {
                    
                    pdb_atom->topo_type = psf_atom.type;
                    pdb_atom->topo_charge = psf_atom.charge;
                    pdb_atom->topo_mass = psf_atom.mass;
                    pdb_atom->chain = psf_atom.segment.empty() ? ' ' : psf_atom.segment[0];
                    updated++;
                    break;
                }
            }
        }
    }

    return updated;
}

std::map<std::string, std::set<std::string>> PSFParser::get_missing_topology_info(
    const std::vector<PDBAtom>& atoms) const {
    std::map<std::string, std::set<std::string>> missing_info;
    
    // First check which residues don't exist in PSF file
    for (const auto& atom : atoms) {
        if (atom_index_.find(atom.residue) == atom_index_.end()) {
            missing_info[atom.residue].insert(atom.name);
            continue;
        }
    }
    
    // Create a copy of the atoms to update with topology info
    std::vector<PDBAtom> atoms_copy = atoms;
    update_pdb_atoms(atoms_copy);
    
    // Check which atoms are still missing topology info
    for (size_t i = 0; i < atoms.size(); ++i) {
        if (!atoms_copy[i].has_topology_info() && 
            atom_index_.find(atoms[i].residue) != atom_index_.end()) {
            missing_info[atoms[i].residue].insert(atoms[i].name);
        }
    }
    
    return missing_info;
}

int PSFParser::update_pdb_atoms_by_order(std::vector<PDBAtom*>& atoms) {
    // Group atoms by residue sequence number
    std::map<std::pair<std::string, int>, std::vector<PDBAtom*>> pdb_atoms_by_res;
    std::map<std::pair<std::string, int>, std::vector<const PSFAtom*>> psf_atoms_by_res;
    
    // Group PDB atoms
    for (auto* atom : atoms) {
        if (atom) {
            pdb_atoms_by_res[std::make_pair(atom->residue, atom->sequence)].push_back(atom);
        }
    }
    
    // Group PSF atoms
    for (const auto& atom : atoms_) {
        psf_atoms_by_res[std::make_pair(atom.residue, atom.residue_number)].push_back(&atom);
    }
    
    // For each residue, map topology data based on order
    int updated_count = 0;
    for (auto& [res_key, pdb_atoms] : pdb_atoms_by_res) {
        auto psf_it = psf_atoms_by_res.find(res_key);
        if (psf_it == psf_atoms_by_res.end()) continue;
        
        auto& psf_atoms = psf_it->second;
        
        // Sort both PDB and PSF atoms by their serial numbers within the residue
        std::sort(pdb_atoms.begin(), pdb_atoms.end(),
            [](const PDBAtom* a, const PDBAtom* b) { return a->serial < b->serial; });
        std::sort(psf_atoms.begin(), psf_atoms.end(),
            [](const PSFAtom* a, const PSFAtom* b) { return a->id < b->id; });
        
        // Map topology data in order
        size_t num_atoms = std::min(pdb_atoms.size(), psf_atoms.size());
        for (size_t i = 0; i < num_atoms; ++i) {
            pdb_atoms[i]->topo_type = psf_atoms[i]->type;
            pdb_atoms[i]->topo_charge = psf_atoms[i]->charge;
            pdb_atoms[i]->topo_mass = psf_atoms[i]->mass;
            pdb_atoms[i]->chain = psf_atoms[i]->segment.empty() ? ' ' : psf_atoms[i]->segment[0];
            updated_count++;
        }
    }
    
    return updated_count;
}

} // namespace io
} // namespace core
} // namespace pygcmc
