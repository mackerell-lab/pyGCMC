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
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
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
            std::cout << "Found " << natoms << " atoms" << std::endl;
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

    return parse_atoms_section(atom_lines);
}

bool PSFParser::parse_atoms_section(const std::vector<std::string>& lines) {
    atoms_.clear();
    atom_index_.clear();

    for (const auto& line : lines) {
        std::istringstream iss(line);
        int atom_id;
        PSFAtom atom;

        // PSF format: ID SEGID RESID RESNAME ATOMNAME ATOMTYPE CHARGE MASS
        if (!(iss >> atom_id >> atom.segment >> atom.residue_number >> atom.residue 
              >> atom.name >> atom.type >> atom.charge >> atom.mass)) {
            std::cerr << "Failed to parse atom line: " << line << std::endl;
            continue;
        }

        std::cout << "Parsed atom: " << atom.residue << " " << atom.name 
                  << " type=" << atom.type 
                  << " charge=" << atom.charge 
                  << " mass=" << atom.mass << std::endl;

        // Add atom to vector and update index
        size_t idx = atoms_.size();
        atoms_.push_back(atom);
        atom_index_[atom.residue][atom.residue_number][atom.name] = idx;
    }

    return !atoms_.empty();
}

bool PSFParser::get_atom_properties(const std::string& residue_name,
                                  const std::string& atom_name,
                                  double& charge,
                                  double& mass) const {
    auto res_it = atom_index_.find(residue_name);
    if (res_it == atom_index_.end()) {
        std::cerr << "Residue not found: " << residue_name << std::endl;
        return false;
    }

    // Get the last defined residue number for this residue type
    auto last_res_num_it = res_it->second.rbegin();
    if (last_res_num_it == res_it->second.rend()) {
        std::cerr << "No residue numbers found for: " << residue_name << std::endl;
        return false;
    }

    auto atom_it = last_res_num_it->second.find(atom_name);
    if (atom_it == last_res_num_it->second.end()) {
        std::cerr << "Atom not found: " << residue_name << " " << atom_name << std::endl;
        return false;
    }

    const PSFAtom& atom = atoms_[atom_it->second];
    charge = atom.charge;
    mass = atom.mass;
    std::cout << "Found atom: " << residue_name << " " << atom_name 
              << " charge=" << charge << " mass=" << mass << std::endl;
    return true;
}

int PSFParser::update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const {
    int updated = 0;
    for (auto& pdb_atom : pdb_atoms) {
        auto res_it = atom_index_.find(pdb_atom.residue);
        if (res_it == atom_index_.end()) {
            std::cerr << "Residue not found: " << pdb_atom.residue << std::endl;
            continue;
        }

        // Find the matching residue number
        // First try exact match
        auto res_num_it = res_it->second.find(pdb_atom.sequence);
        if (res_num_it == res_it->second.end()) {
            // If exact match fails, try to find any matching residue number
            // that has the same atom name
            bool found = false;
            for (const auto& [res_num, atoms] : res_it->second) {
                auto atom_it = atoms.find(pdb_atom.name);
                if (atom_it != atoms.end()) {
                    res_num_it = res_it->second.find(res_num);
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Residue number not found: " << pdb_atom.residue << " " << pdb_atom.sequence << std::endl;
                continue;
            }
        }

        auto atom_it = res_num_it->second.find(pdb_atom.name);
        if (atom_it == res_num_it->second.end()) {
            std::cerr << "Atom not found: " << pdb_atom.residue << " " << pdb_atom.name << std::endl;
            continue;
        }

        const PSFAtom& atom = atoms_[atom_it->second];
        pdb_atom.topo_type = atom.type;
        pdb_atom.topo_charge = atom.charge;
        pdb_atom.topo_mass = atom.mass;
        std::cout << "Updated atom: " << pdb_atom.residue << " " << pdb_atom.name 
                  << " type=" << pdb_atom.topo_type 
                  << " charge=" << pdb_atom.topo_charge 
                  << " mass=" << pdb_atom.topo_mass << std::endl;
        updated++;
    }
    return updated;
}

std::map<std::string, std::set<std::string>> PSFParser::get_missing_topology_info(
    const std::vector<PDBAtom>& atoms) const {
    std::map<std::string, std::set<std::string>> missing_info;
    
    // Create a copy of the atoms to update with topology info
    std::vector<PDBAtom> atoms_copy = atoms;
    update_pdb_atoms(atoms_copy);
    
    // Check which atoms are still missing topology info
    for (size_t i = 0; i < atoms.size(); ++i) {
        if (!atoms_copy[i].has_topology_info()) {
            missing_info[atoms[i].residue].insert(atoms[i].name);
        }
    }
    
    return missing_info;
}

} // namespace io
} // namespace core
} // namespace pygcmc
