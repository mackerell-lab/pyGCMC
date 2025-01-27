#include "psfParser.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace io {

bool PSFParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open PSF file: " << filename << std::endl;
        return false;
    }

    // Read and validate PSF header
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Failed to read PSF header" << std::endl;
        return false;
    }

    // Check for PSF header
    if (line.find("PSF") == std::string::npos) {
        std::cerr << "Invalid PSF file: missing PSF header" << std::endl;
        return false;
    }

    // Parse title section
    if (!parse_title(file, topology)) {
        std::cerr << "Failed to parse title section" << std::endl;
        return false;
    }

    // Parse atoms section
    if (!parse_atoms(file, topology)) {
        std::cerr << "Failed to parse atoms section" << std::endl;
        return false;
    }

    return true;
}

bool PSFParser::read_section_header(std::ifstream& file, const std::string& expected_header, int& count) {
    std::string line;
    std::getline(file, line);
    if (line.empty()) return false;
    
    // Verify the header if provided
    if (!expected_header.empty() && line.find(expected_header) == std::string::npos) {
        return false;
    }
    
    std::istringstream iss(line);
    iss >> count;
    return true;
}

std::vector<int> PSFParser::read_index_block(std::ifstream& file, int expected_count, int indices_per_item) {
    std::vector<int> indices;
    indices.reserve(static_cast<size_t>(expected_count) * static_cast<size_t>(indices_per_item));
    
    std::string line;
    while (indices.size() < static_cast<size_t>(expected_count) * static_cast<size_t>(indices_per_item)) {
        std::getline(file, line);
        std::istringstream iss(line);
        int idx;
        while (iss >> idx) {
            indices.push_back(idx - 1);  // Convert to 0-based indexing
        }
    }
    return indices;
}

bool PSFParser::parse_title(std::ifstream& file, model::Topology& topology) {
    std::string line;
    int num_titles = 0;

    // Read lines until we find !NTITLE
    while (std::getline(file, line)) {
        std::string trimmed = trim(line);
        if (trimmed.empty()) {
            continue;
        }
        if (trimmed.find("!NTITLE") != std::string::npos) {
            std::istringstream iss(trimmed);
            iss >> num_titles;
            break;
        }
    }

    // Read exactly num_titles lines
    for (int i = 0; i < num_titles; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Failed to read title lines" << std::endl;
            return false;
        }
        if (!line.empty()) {
            topology.add_title(line);
        }
    }
    return true;
}

bool PSFParser::parse_atoms(std::ifstream& file, model::Topology& topology) {
    std::string line;
    int num_atoms = 0;

    // Step 1: Find the line containing "!NATOM" and parse the integer before it
    while (std::getline(file, line)) {
        std::string trimmedLine = trim(line);
        if (trimmedLine.empty()) {
            continue;
        }
        if (trimmedLine.find("!NATOM") != std::string::npos) {
            std::istringstream iss(trimmedLine);
            iss >> num_atoms;
            break;
        }
    }

    if (num_atoms <= 0) {
        std::cerr << "Invalid number of atoms: " << num_atoms << std::endl;
        return false;
    }

    // Step 2: Read exactly num_atoms lines
    for (int i = 0; i < num_atoms; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Failed to read atom line " << i + 1 << std::endl;
            return false;
        }

        // Split by whitespace
        std::istringstream iss(line);

        // Typical PSF has these 9 fields in each ATOM line (XPLOR style):
        //  1) atomIndex (int)
        //  2) segmentName (string)
        //  3) residueNumber (int)
        //  4) residueName (string)
        //  5) atomName (string)
        //  6) atomType (string)
        //  7) charge (double)
        //  8) mass (double)
        //  9) extra integer (often 0)
        
        int atomIndex;
        std::string segment_name;
        int residue_number;
        std::string residue_name;
        std::string atom_name;
        std::string atom_type;
        double charge;
        double mass;
        int unusedField; // often 0 or some integer

        // Try reading all fields - if some files omit the last integer, we'll handle that gracefully
        bool parsedOK = false;

        // Try reading 9 fields first
        if (iss >> atomIndex >> segment_name >> residue_number
                >> residue_name >> atom_name >> atom_type
                >> charge >> mass >> unusedField) {
            parsedOK = true;
        } else {
            // Clear stream state and rewind
            iss.clear();
            iss.seekg(0);

            // Try reading only 8 fields if the file doesn't have the extra integer
            if (iss >> atomIndex >> segment_name >> residue_number
                    >> residue_name >> atom_name >> atom_type
                    >> charge >> mass) {
                parsedOK = true;
            }
        }

        if (!parsedOK) {
            std::cerr << "Failed to parse atom line " << (i + 1)
                      << ": " << line << std::endl;
            return false;
        }

        // Add to topology (ignore atomIndex and unusedField)
        topology.add_atom(
            atom_name,
            atom_type,
            charge,
            mass,
            residue_name,
            residue_number,
            segment_name
        );
    }

    return true;
}

// Helper function to trim whitespace
std::string PSFParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

bool PSFParser::parse_bonds(std::ifstream& file, model::Topology& topology) {
    std::string line;
    std::getline(file, line);  // Get the number of bonds
    int num_bonds = std::stoi(line);

    std::vector<int> bond_indices;
    bond_indices.reserve(static_cast<size_t>(num_bonds) * 2);
    while (bond_indices.size() < static_cast<size_t>(num_bonds) * 2) {
        std::getline(file, line);
        std::istringstream iss(line);
        int idx;
        while (iss >> idx) {
            bond_indices.push_back(idx - 1);  // Convert to 0-based indexing
        }
    }

    for (size_t i = 0; i < bond_indices.size(); i += 2) {
        topology.add_bond(bond_indices[i], bond_indices[i + 1]);
    }
    return true;
}

bool PSFParser::parse_angles(std::ifstream& file, model::Topology& topology) {
    std::string line;
    std::getline(file, line);  // Get the number of angles
    int num_angles = std::stoi(line);

    std::vector<int> angle_indices;
    angle_indices.reserve(static_cast<size_t>(num_angles) * 3);
    while (angle_indices.size() < static_cast<size_t>(num_angles) * 3) {
        std::getline(file, line);
        std::istringstream iss(line);
        int idx;
        while (iss >> idx) {
            angle_indices.push_back(idx - 1);  // Convert to 0-based indexing
        }
    }

    for (size_t i = 0; i < angle_indices.size(); i += 3) {
        topology.add_angle(angle_indices[i], angle_indices[i + 1], angle_indices[i + 2]);
    }
    return true;
}

bool PSFParser::parse_dihedrals(std::ifstream& file, model::Topology& topology) {
    std::string line;
    std::getline(file, line);  // Get the number of dihedrals
    int num_dihedrals = std::stoi(line);

    std::vector<int> dihedral_indices;
    dihedral_indices.reserve(static_cast<size_t>(num_dihedrals) * 4);
    while (dihedral_indices.size() < static_cast<size_t>(num_dihedrals) * 4) {
        std::getline(file, line);
        std::istringstream iss(line);
        int idx;
        while (iss >> idx) {
            dihedral_indices.push_back(idx - 1);  // Convert to 0-based indexing
        }
    }

    for (size_t i = 0; i < dihedral_indices.size(); i += 4) {
        topology.add_dihedral(
            dihedral_indices[i],
            dihedral_indices[i + 1],
            dihedral_indices[i + 2],
            dihedral_indices[i + 3]
        );
    }
    return true;
}

bool PSFParser::parse_impropers(std::ifstream& file, model::Topology& topology) {
    int num_impropers;
    if (!read_section_header(file, "!NIMPHI", num_impropers)) return false;

    auto indices = read_index_block(file, num_impropers, 4);
    for (size_t i = 0; i < indices.size(); i += 4) {
        topology.add_improper(indices[i], indices[i + 1], indices[i + 2], indices[i + 3]);
    }
    return true;
}

bool PSFParser::parse_donors(std::ifstream& file, model::Topology& topology) {
    int num_donors;
    if (!read_section_header(file, "!NDON", num_donors)) return false;

    auto indices = read_index_block(file, num_donors, 2);
    for (size_t i = 0; i < indices.size(); i += 2) {
        topology.add_donor(indices[i], indices[i + 1]);
    }
    return true;
}

bool PSFParser::parse_acceptors(std::ifstream& file, model::Topology& topology) {
    int num_acceptors;
    if (!read_section_header(file, "!NACC", num_acceptors)) return false;

    auto indices = read_index_block(file, num_acceptors, 1);
    for (const auto& idx : indices) {
        topology.add_acceptor(idx);
    }
    return true;
}

bool PSFParser::parse_nonbonded_exclusions(std::ifstream& file, model::Topology& topology) {
    int num_atoms;
    if (!read_section_header(file, "!NNB", num_atoms)) return false;

    std::string line;
    for (int i = 0; i < num_atoms; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        int atom_idx, num_exclusions;
        iss >> atom_idx >> num_exclusions;
        atom_idx--;  // Convert to 0-based indexing
        
        for (int j = 0; j < num_exclusions; ++j) {
            int excluded_atom;
            iss >> excluded_atom;
            topology.add_nonbonded_exclusion(atom_idx, excluded_atom - 1);
        }
    }
    return true;
}

bool PSFParser::parse_groups(std::ifstream& file, model::Topology& topology) {
    int num_groups;
    if (!read_section_header(file, "!NGRP", num_groups)) return false;

    std::string line;
    for (int i = 0; i < num_groups; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        int group_id;
        std::vector<int> atoms;
        std::string type;
        
        iss >> group_id;
        int atom_idx;
        while (iss >> atom_idx) {
            atoms.push_back(atom_idx - 1);
        }
        
        topology.add_group(group_id, atoms, type);
    }
    return true;
}

bool PSFParser::parse_cmap(std::ifstream& file, model::Topology& topology) {
    int num_cmap;
    if (!read_section_header(file, "!NCRTERM", num_cmap)) return false;

    auto indices = read_index_block(file, num_cmap, 8);
    for (size_t i = 0; i < indices.size(); i += 8) {
        std::array<int, 8> cmap_atoms;
        std::copy(indices.begin() + i, indices.begin() + i + 8, cmap_atoms.begin());
        topology.add_cmap(cmap_atoms);
    }
    return true;
}

} // namespace io
} // namespace pygcmc
