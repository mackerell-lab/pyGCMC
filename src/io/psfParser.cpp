#include "psfParser.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <array>

namespace pygcmc {
namespace io {

bool PSFParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    std::getline(file, line);  // Read "PSF" header
    if (line.find("PSF") == std::string::npos) {
        return false;  // Not a PSF file
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '*') continue;

        if (line.find("!NTITLE") != std::string::npos) {
            if (!parse_title(file, topology)) return false;
        }
        else if (line.find("!NATOM") != std::string::npos) {
            if (!parse_atoms(file, topology)) return false;
        }
        else if (line.find("!NBOND") != std::string::npos) {
            if (!parse_bonds(file, topology)) return false;
        }
        else if (line.find("!NTHETA") != std::string::npos) {
            if (!parse_angles(file, topology)) return false;
        }
        else if (line.find("!NPHI") != std::string::npos) {
            if (!parse_dihedrals(file, topology)) return false;
        }
        else if (line.find("!NIMPHI") != std::string::npos) {
            if (!parse_impropers(file, topology)) return false;
        }
        else if (line.find("!NDON") != std::string::npos) {
            if (!parse_donors(file, topology)) return false;
        }
        else if (line.find("!NACC") != std::string::npos) {
            if (!parse_acceptors(file, topology)) return false;
        }
        else if (line.find("!NNB") != std::string::npos) {
            if (!parse_nonbonded_exclusions(file, topology)) return false;
        }
        else if (line.find("!NGRP") != std::string::npos) {
            if (!parse_groups(file, topology)) return false;
        }
        else if (line.find("!NCRTERM") != std::string::npos) {
            if (!parse_cmap(file, topology)) return false;
        }
    }

    return true;
}

bool PSFParser::read_section_header(std::ifstream& file, const std::string& expected_header, int& count) {
    std::string line;
    std::getline(file, line);
    if (line.empty()) return false;
    
    std::istringstream iss(line);
    iss >> count;
    return true;
}

std::vector<int> PSFParser::read_index_block(std::ifstream& file, int expected_count, int indices_per_item) {
    std::vector<int> indices;
    indices.reserve(expected_count * indices_per_item);
    
    std::string line;
    while (indices.size() < expected_count * indices_per_item) {
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
    int num_titles;
    if (!read_section_header(file, "!NTITLE", num_titles)) return false;

    std::string line;
    for (int i = 0; i < num_titles; ++i) {
        std::getline(file, line);
        if (!line.empty()) {
            topology.add_title(line);
        }
    }
    return true;
}

bool PSFParser::parse_atoms(std::ifstream& file, model::Topology& topology) {
    std::string line;
    std::getline(file, line);  // Get the number of atoms
    int num_atoms = std::stoi(line);

    for (int i = 0; i < num_atoms; ++i) {
        std::getline(file, line);
        if (line.empty()) continue;

        // PSF atom format:
        // atomid segname resid resname atomname atomtype charge mass
        std::istringstream iss(line);
        int atomid;
        std::string segname, resname, atomname, atomtype;
        int resid;
        double charge, mass;

        iss >> atomid >> segname >> resid >> resname >> atomname >> atomtype >> charge >> mass;

        topology.add_atom(atomname, atomtype, charge, mass, resname, resid, segname);
    }
    return true;
}

bool PSFParser::parse_bonds(std::ifstream& file, model::Topology& topology) {
    std::string line;
    std::getline(file, line);  // Get the number of bonds
    int num_bonds = std::stoi(line);

    std::vector<int> bond_indices;
    while (bond_indices.size() < num_bonds * 2) {
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
    while (angle_indices.size() < num_angles * 3) {
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
    while (dihedral_indices.size() < num_dihedrals * 4) {
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
    for (int idx : indices) {
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
