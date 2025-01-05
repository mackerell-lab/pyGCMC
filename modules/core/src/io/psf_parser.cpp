// modules/core/src/io/psf_parser.cpp

#include "pygcmc/core/io/psf_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <unordered_map>

namespace pygcmc {
namespace core {
namespace io {

bool PSFParser::parse(const std::string& filename, PSFParsingMode mode) {
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
                if (mode == PSFParsingMode::Exact) {
                    throw std::runtime_error("Invalid number of atoms in PSF file: " + std::to_string(natoms));
                } else {
                    // In rough mode, proceed without atom count
                    natoms = 0;
                }
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
        if (mode == PSFParsingMode::Exact) {
            throw std::runtime_error("No atoms found in PSF file");
        } else {
            // In rough mode, it's acceptable to have no atoms
            return true;
        }
    }

    if (mode == PSFParsingMode::Exact) {
        return parse_atoms_section(atom_lines, mode);
    } else {
        return parse_atoms_section_rough(atom_lines);
    }
}

bool PSFParser::parse_atoms_section(const std::vector<std::string>& lines, PSFParsingMode mode) {
    // Only clear existing atoms if this is the first file being parsed
    if (is_first_file_ && mode == PSFParsingMode::Exact) {
        atoms_.clear();
        atom_index_.clear();
        id_index_.clear();
        is_first_file_ = false;
    }

    for (const auto& line : lines) {
        std::istringstream iss(line);
        PSFAtom atom;

        // PSF exact format: ID SEGID RESID RESNAME ATOMNAME ATOMTYPE CHARGE MASS
        if (!(iss >> atom.id 
                  >> atom.segment 
                  >> atom.residue_number 
                  >> atom.residue 
                  >> atom.name 
                  >> atom.type 
                  >> atom.charge 
                  >> atom.mass)) 
        {
            if (mode == PSFParsingMode::Exact) {
                std::cerr << "Failed to parse atom line in exact mode: " << line << std::endl;
                throw std::runtime_error("Exact parsing failed due to malformed atom line.");
            } else {
                std::cerr << "Failed to parse atom line in exact mode: " << line << std::endl;
                continue; // Skip malformed lines in exact mode
            }
        }

        // Store the atom
        size_t idx = atoms_.size();
        atoms_.push_back(atom);

        // Update indexes
        atom_index_[atom.residue][atom.residue_number][atom.name] = idx;
        id_index_[atom.id] = idx;
    }

    return !atoms_.empty();
}

bool PSFParser::parse_atoms_section_rough(const std::vector<std::string>& lines) {
    for (const auto& line : lines) {
        std::istringstream iss(line);
        PSFAtom atom;
        std::string dummy; // For skipping fields we don't need

        // Try to extract the essential fields
        if (!(iss >> dummy           // Skip ID
                  >> atom.segment    // Keep segment
                  >> dummy           // Skip residue number
                  >> atom.residue    // Keep residue name
                  >> atom.name       // Keep atom name
                  >> atom.type))     // Keep atom type
        {
            std::cerr << "Failed to parse atom line in rough mode: " << line << std::endl;
            continue; // Skip malformed lines
        }

        // Assign default values for missing fields
        atom.id = 0;  // Not reliable in rough mode
        atom.residue_number = 0;  // Not reliable in rough mode
        atom.charge = 0.0;  // Will be assigned later if needed
        atom.mass = 0.0;   // Will be assigned later if needed

        // Store the atom
        size_t idx = atoms_.size();
        atoms_.push_back(atom);

        // In rough mode, we only index by residue name and atom name
        atom_index_[atom.residue][0][atom.name] = idx;
    }

    return !atoms_.empty();
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

    // Try exact residue number first
    auto res_num_it = res_it->second.find(residue_number);
    if (res_num_it != res_it->second.end()) {
        auto atom_it = res_num_it->second.find(atom_name);
        if (atom_it != res_num_it->second.end()) {
            const PSFAtom& atom = atoms_[atom_it->second];
            charge = atom.charge;
            mass = atom.mass;
            return true;
        }
    }

    // If not found with exact residue number and we're not looking for a specific number,
    // try to find any matching atom
    if (residue_number == 0) {
        for (const auto& [num, atoms] : res_it->second) {
            auto atom_it = atoms.find(atom_name);
            if (atom_it != atoms.end()) {
                const PSFAtom& atom = atoms_[atom_it->second];
                charge = atom.charge;
                mass = atom.mass;
                return true;
            }
        }
    }

    std::cerr << "[get_atom_properties] Atom not found: "
              << residue_name << " " << residue_number << " " << atom_name << std::endl;
    return false;
}

bool PSFParser::get_atom_properties(const std::string& residue_name,
                                  const std::string& atom_name,
                                  double& charge,
                                  double& mass) const {
    return get_atom_properties(residue_name, 0, atom_name, charge, mass);
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
    
    // Create a map for faster lookup of PSF atoms
    std::unordered_map<std::string, std::unordered_map<int, std::unordered_map<std::string, const PSFAtom*>>> psf_atom_map;
    for (const auto& psf_atom : atoms_) {
        psf_atom_map[psf_atom.residue][psf_atom.residue_number][psf_atom.name] = &psf_atom;
    }
    
    // First pass: match by exact residue number and atom name
    for (auto* pdb_atom : pdb_atoms) {
        if (!pdb_atom) continue;
        
        // Try exact match first
        auto res_it = psf_atom_map.find(pdb_atom->residue);
        if (res_it != psf_atom_map.end()) {
            auto seq_it = res_it->second.find(pdb_atom->sequence);
            if (seq_it != res_it->second.end()) {
                auto atom_it = seq_it->second.find(pdb_atom->name);
                if (atom_it != seq_it->second.end()) {
                    const PSFAtom* psf_atom = atom_it->second;
                    // Only update if not already set or if we have a better match
                    if (pdb_atom->topo_type.empty() || 
                        (psf_atom->residue_number == pdb_atom->sequence && 
                         psf_atom->name == pdb_atom->name)) {
                        pdb_atom->topo_type = psf_atom->type;
                        pdb_atom->topo_charge = psf_atom->charge;
                        pdb_atom->topo_mass = psf_atom->mass;
                        pdb_atom->chain = psf_atom->segment.empty() ? ' ' : psf_atom->segment[0];
                        updated++;
                    }
                }
            }
        }
    }
    
    // Second pass: match by residue name and atom name only if not matched in first pass
    for (auto* pdb_atom : pdb_atoms) {
        if (!pdb_atom || !pdb_atom->topo_type.empty()) continue;  // Skip if already matched
        
        auto res_it = psf_atom_map.find(pdb_atom->residue);
        if (res_it != psf_atom_map.end()) {
            // Try to find a match with any sequence number
            for (const auto& [seq_num, atom_data] : res_it->second) {
                auto atom_it = atom_data.find(pdb_atom->name);
                if (atom_it != atom_data.end()) {
                    const PSFAtom* psf_atom = atom_it->second;
                    pdb_atom->topo_type = psf_atom->type;
                    pdb_atom->topo_charge = psf_atom->charge;
                    pdb_atom->topo_mass = psf_atom->mass;
                    pdb_atom->chain = psf_atom->segment.empty() ? ' ' : psf_atom->segment[0];
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

int PSFParser::update_pdb_atoms_from_multiple_psf(std::vector<PDBAtom*>& pdb_atoms,
                                                const std::vector<std::string>& psf_files) {
    int total_updated = 0;
    std::vector<PDBAtom*> remaining_atoms = pdb_atoms;

    // Try to update atoms from each PSF file
    for (const auto& psf_file : psf_files) {
        PSFParser parser;
        if (!parser.parse(psf_file)) {
            continue;
        }

        // Try to update all atoms with this PSF file
        int updated = parser.update_pdb_atoms(pdb_atoms);
        if (updated > 0) {
            total_updated += updated;
        }
    }

    return total_updated;
}

int PSFParser::update_pdb_atoms_multi_residue(std::vector<PDBAtom*>& pdb_atoms, const std::string& psf_file) {
    PSFParser parser;
    if (!parser.parse(psf_file, PSFParsingMode::Exact)) {
        return 0;
    }
    return parser.update_pdb_atoms(pdb_atoms);
}

int PSFParser::update_pdb_atoms_single_residue(std::vector<PDBAtom*>& pdb_atoms, const std::string& psf_file, const std::string& target_residue) {
    PSFParser parser;
    if (!parser.parse(psf_file, PSFParsingMode::Exact)) {
        return 0;
    }
    
    int updated = 0;
    for (auto* pdb_atom : pdb_atoms) {
        if (!pdb_atom) continue;
        if (pdb_atom->residue != target_residue) continue;
        
        // Try to find matching atom in PSF file
        for (const auto& psf_atom : parser.atoms_) {
            if (psf_atom.residue == target_residue && psf_atom.name == pdb_atom->name) {
                pdb_atom->topo_type = psf_atom.type;
                pdb_atom->topo_charge = psf_atom.charge;
                pdb_atom->topo_mass = psf_atom.mass;
                pdb_atom->chain = psf_atom.segment.empty() ? ' ' : psf_atom.segment[0];
                updated++;
                break;
            }
        }
    }
    
    return updated;
}

} // namespace io
} // namespace core
} // namespace pygcmc
