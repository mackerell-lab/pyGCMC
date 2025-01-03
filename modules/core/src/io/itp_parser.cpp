// modules/core/src/io/itp_parser.cpp

#include "pygcmc/core/io/itp_parser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace pygcmc {
namespace core {
namespace io {

bool ITPParser::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }

    std::vector<std::string> current_section_lines;
    std::string current_section;
    std::string line;

    while (std::getline(file, line)) {
        // Remove leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        // Skip empty lines and comments
        if (line.empty() || line[0] == ';' || line[0] == '#') {
            continue;
        }

        // Check for section header
        if (line[0] == '[' && line.back() == ']') {
            // Process previous section
            if (!current_section.empty()) {
                if (current_section == "atoms") {
                    if (!parse_atoms_section(current_section_lines)) {
                        std::cerr << "Failed to parse [ atoms ] section" << std::endl;
                        return false;
                    }
                }
            }

            // Set current section
            current_section = line.substr(1, line.size() - 2);
            // Remove leading/trailing whitespace from section name
            current_section.erase(0, current_section.find_first_not_of(" \t"));
            current_section.erase(current_section.find_last_not_of(" \t") + 1);
            current_section_lines.clear();
            continue;
        }

        // Add line to current section
        if (!current_section.empty()) {
            current_section_lines.push_back(line);
        }
    }

    // Process last section
    if (!current_section.empty() && current_section == "atoms") {
        if (!parse_atoms_section(current_section_lines)) {
            std::cerr << "Failed to parse [ atoms ] section" << std::endl;
            return false;
        }
    }

    return true;
}

bool ITPParser::parse_atoms_section(const std::vector<std::string>& lines) {
    if (lines.empty()) {
        std::cerr << "[ atoms ] section is empty" << std::endl;
        return false;
    }

    // Skip header line if present
    size_t start = 0;
    if (lines[0].find("nr") != std::string::npos) {
        start = 1;
    }

    for (size_t i = start; i < lines.size(); ++i) {
        const auto& line = lines[i];
        if (line.empty() || line[0] == ';') {
            continue;
        }

        std::istringstream iss(line);
        int nr, resnr, cgnr;
        std::string type, resname, atom;
        double charge, mass;

        // Parse atom line
        if (!(iss >> nr >> type >> resnr >> resname >> atom >> cgnr >> charge >> mass)) {
            std::cerr << "Failed to parse atom line: " << line << std::endl;
            continue;
        }

        ITPAtom itp_atom;
        itp_atom.name = atom;
        itp_atom.type = type;
        itp_atom.resname = resname;
        itp_atom.resid = resnr;
        itp_atom.charge = charge;
        itp_atom.mass = mass;

        itp_atoms_.push_back(itp_atom);
        ResidueKey key{resname, resnr};
        atom_index_[key][atom] = itp_atoms_.size() - 1;
    }

    return !itp_atoms_.empty();
}

bool ITPParser::get_atom_properties(const std::string& residue_name,
                                  const std::string& atom_name,
                                  double& charge,
                                  double& mass) const {
    // Try residue number 1 first (most common case)
    ResidueKey key{residue_name, 1};
    auto res_it = atom_index_.find(key);
    if (res_it == atom_index_.end()) {
        // If not found, try other residue numbers
        for (const auto& [key, atoms] : atom_index_) {
            if (key.resname == residue_name) {
                auto atom_it = atoms.find(atom_name);
                if (atom_it != atoms.end()) {
                    const ITPAtom& atom = itp_atoms_[atom_it->second];
                    charge = atom.charge;
                    mass = atom.mass;
                    return true;
                }
            }
        }
        return false;
    }

    auto atom_it = res_it->second.find(atom_name);
    if (atom_it == res_it->second.end()) {
        return false;
    }

    const ITPAtom& atom = itp_atoms_[atom_it->second];
    charge = atom.charge;
    mass = atom.mass;
    return true;
}

int ITPParser::update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const {
    int updated = 0;
    
    for (auto& pdb_atom : pdb_atoms) {
        // First try to find the atom in the first residue (usually residue 1)
        ResidueKey key{pdb_atom.residue, 1};
        auto res_it = atom_index_.find(key);
        
        if (res_it == atom_index_.end()) {
            continue;
        }

        auto atom_it = res_it->second.find(pdb_atom.name);
        if (atom_it == res_it->second.end()) {
            continue;
        }

        const ITPAtom& itp_atom = itp_atoms_[atom_it->second];
        pdb_atom.topo_type = itp_atom.type;
        pdb_atom.topo_charge = itp_atom.charge;
        pdb_atom.topo_mass = itp_atom.mass;
        updated++;
    }
    return updated;
}

int ITPParser::update_pdb_atoms(std::vector<PDBAtom*>& pdb_atoms) const {
    int updated = 0;
    
    for (auto* pdb_atom : pdb_atoms) {
        if (!pdb_atom) continue;
        
        // First try to find the atom in the first residue (usually residue 1)
        ResidueKey key{pdb_atom->residue, 1};
        auto res_it = atom_index_.find(key);
        
        if (res_it == atom_index_.end()) {
            continue;
        }

        auto atom_it = res_it->second.find(pdb_atom->name);
        if (atom_it == res_it->second.end()) {
            continue;
        }

        const ITPAtom& itp_atom = itp_atoms_[atom_it->second];
        pdb_atom->topo_type = itp_atom.type;
        pdb_atom->topo_charge = itp_atom.charge;
        pdb_atom->topo_mass = itp_atom.mass;
        updated++;
    }
    return updated;
}

std::map<std::string, std::set<std::string>> ITPParser::get_missing_topology_info(
    const std::vector<PDBAtom>& atoms) const {
    std::map<std::string, std::set<std::string>> missing_info;
    
    for (const auto& atom : atoms) {
        ResidueKey key{atom.residue, 1};
        auto res_it = atom_index_.find(key);
        if (res_it == atom_index_.end()) {
            missing_info[atom.residue].insert(atom.name);
            continue;
        }

        auto atom_it = res_it->second.find(atom.name);
        if (atom_it == res_it->second.end()) {
            missing_info[atom.residue].insert(atom.name);
        }
    }
    
    return missing_info;
}

} // namespace io
} // namespace core
} // namespace pygcmc

