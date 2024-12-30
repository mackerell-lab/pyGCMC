// modules/core/src/io/top_parser.cpp

#include "pygcmc/core/io/top_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace pygcmc {
namespace core {
namespace io {

bool TopParser::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }

    std::vector<std::string> atom_lines;
    bool in_atoms_section = false;
    std::string line;

    while (std::getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == ';') {
            continue;
        }

        // Remove leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (line == "[ atoms ]") {
            std::cout << "Found atoms section" << std::endl;
            in_atoms_section = true;
            continue;
        } else if (line[0] == '[' && in_atoms_section) {
            // End of atoms section
            break;
        }

        if (in_atoms_section && !line.empty()) {
            atom_lines.push_back(line);
        }
    }

    std::cout << "Found " << atom_lines.size() << " atom lines" << std::endl;
    return parse_atoms_section(atom_lines);
}

bool TopParser::parse_atoms_section(const std::vector<std::string>& lines) {
    atoms_.clear();
    atom_index_.clear();

    for (const auto& line : lines) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;

        // Split line into tokens
        while (iss >> token) {
            tokens.push_back(token);
        }

        // Skip lines that don't have enough tokens
        // Format: nr type resnr residue atom cgnr charge mass
        if (tokens.size() < 8) {
            continue;
        }

        TopAtom atom;
        atom.type = tokens[1];      // Atom type from topology
        try {
            atom.residue_number = std::stoi(tokens[2]);  // Parse residue number
        } catch (const std::exception&) {
            std::cerr << "Failed to parse residue number for line: " << line << std::endl;
            continue;
        }
        atom.residue = tokens[3];
        atom.name = tokens[4];
        
        try {
            atom.charge = std::stod(tokens[6]);
            atom.mass = std::stod(tokens[7]);
            std::cout << "Parsed atom: " << atom.residue << " " << atom.name 
                      << " type=" << atom.type 
                      << " charge=" << atom.charge 
                      << " mass=" << atom.mass << std::endl;
        } catch (const std::exception&) {
            std::cerr << "Failed to parse charge/mass for line: " << line << std::endl;
            continue;
        }

        // Add atom to vector and update index
        size_t idx = atoms_.size();
        atoms_.push_back(atom);
        // Use residue number in the index
        atom_index_[atom.residue][atom.residue_number][atom.name] = idx;
    }

    return !atoms_.empty();
}

bool TopParser::get_atom_properties(const std::string& residue_name,
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

    const TopAtom& atom = atoms_[atom_it->second];
    charge = atom.charge;
    mass = atom.mass;
    std::cout << "Found atom: " << residue_name << " " << atom_name 
              << " charge=" << charge << " mass=" << mass << std::endl;
    return true;
}

int TopParser::update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const {
    int updated = 0;
    for (auto& pdb_atom : pdb_atoms) {
        auto res_it = atom_index_.find(pdb_atom.residue);
        if (res_it == atom_index_.end()) {
            std::cerr << "Residue not found: " << pdb_atom.residue << std::endl;
            continue;
        }

        // Find the matching residue number
        auto res_num_it = res_it->second.find(pdb_atom.sequence);
        if (res_num_it == res_it->second.end()) {
            std::cerr << "Residue number not found: " << pdb_atom.residue << " " << pdb_atom.sequence << std::endl;
            continue;
        }

        auto atom_it = res_num_it->second.find(pdb_atom.name);
        if (atom_it == res_num_it->second.end()) {
            std::cerr << "Atom not found: " << pdb_atom.residue << " " << pdb_atom.name << std::endl;
            continue;
        }

        const TopAtom& atom = atoms_[atom_it->second];
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

std::map<std::string, std::set<std::string>> TopParser::get_missing_topology_info(
    const std::vector<PDBAtom>& atoms) const {
    std::map<std::string, std::set<std::string>> missing_info;
    
    for (const auto& atom : atoms) {
        if (!atom.has_topology_info()) {
            missing_info[atom.residue].insert(atom.name);
        }
    }
    
    return missing_info;
}

} // namespace io
} // namespace core
} // namespace pygcmc
