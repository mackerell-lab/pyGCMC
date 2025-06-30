// src/io/topology/psfParserSectionsBasic.cpp

#include "psfParserSectionsBasic.hpp"
#include "psfParserStringUtils.hpp"
#include <sstream>
#include <iostream>
#include <algorithm>

namespace pygcmc {
namespace io {

bool PSFParserSectionsBasic::parse_title_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology) {
    // Parse number of titles from the current line
    std::istringstream iss(PSFParserStringUtils::trim(lines[current_line]));
    int num_titles = 0;
    if (!(iss >> num_titles)) {
        std::cerr << "Failed to parse number of title lines" << std::endl;
        return false;
    }

    // Read exactly num_titles lines
    for (int i = 0; i < num_titles; ++i) {
        current_line++;
        if (current_line >= lines.size()) {
            std::cerr << "Failed to read title line " << i + 1 << std::endl;
            return false;
        }
        if (!lines[current_line].empty()) {
            topology.add_title(lines[current_line]);
        }
    }
    return true;
}

bool PSFParserSectionsBasic::parse_atoms_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    // Parse number of atoms from the first line (header)
    std::istringstream iss(PSFParserStringUtils::trim(lines[0]));  // Always use first line as header
    int num_atoms = 0;
    std::string marker;  // For "!NATOM" marker
    
    // Try to parse the line with or without the marker
    if (!(iss >> num_atoms)) {
        iss.clear();
        iss.seekg(0);
        if (!(iss >> num_atoms >> marker)) {
            std::cerr << "Failed to parse number of atoms" << std::endl;
            return false;
        }
    }

    if (num_atoms <= 0) {
        std::cerr << "Invalid number of atoms: " << num_atoms << std::endl;
        return false;
    }

    // Start from line 1 (after header) and read atom lines
    int atoms_read = 0;
    for (size_t i = 1; i < lines.size() && atoms_read < num_atoms; ++i) {
        std::string line = PSFParserStringUtils::trim(lines[i]);
        if (line.empty()) continue;

        std::istringstream iss(line);
        int atomIndex;
        std::string segment_name;
        int residue_number;
        std::string residue_name;
        std::string atom_name;
        std::string atom_type;
        double charge;
        double mass;
        int unusedField = 0;  // Optional field

        // Try to parse with all fields
        if (iss >> atomIndex >> segment_name >> residue_number
                >> residue_name >> atom_name >> atom_type
                >> charge >> mass >> unusedField) {
            // Successfully parsed all fields
        }
        // Try without the unused field
        else {
            iss.clear();
            iss.seekg(0);
            if (!(iss >> atomIndex >> segment_name >> residue_number
                    >> residue_name >> atom_name >> atom_type
                    >> charge >> mass)) {
                std::cerr << "Failed to parse atom line: " << line << std::endl;
                return false;
            }
        }

        // Convert to 0-based indexing
        atomIndex--;
        if (atomIndex != atoms_read) {
            std::cerr << "Atom index mismatch at line " << (atoms_read + 1) << std::endl;
            return false;
        }

        topology.add_atom(
            atom_name,
            atom_type,
            charge,
            mass,
            residue_name,
            residue_number,
            segment_name
        );
        atoms_read++;
    }

    if (atoms_read != num_atoms) {
        std::cerr << "Expected " << num_atoms << " atoms but read " << atoms_read << std::endl;
        return false;
    }

    return true;
}

bool PSFParserSectionsBasic::parse_bonds_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    // Parse number of bonds from the first line
    std::istringstream iss(PSFParserStringUtils::trim(lines[0]));
    int num_bonds = 0;
    std::string marker;  // For "!NBOND" marker
    
    // Try to parse the line with or without the marker
    if (!(iss >> num_bonds)) {
        iss.clear();
        iss.seekg(0);
        if (!(iss >> num_bonds >> marker)) {
            std::cerr << "Failed to parse number of bonds" << std::endl;
            return false;
        }
    }

    if (num_bonds <= 0) {
        // Bonds are optional, so just return true if there are none
        return true;
    }

    // Read bond indices
    std::vector<int> bond_indices;
    bond_indices.reserve(num_bonds * 2);

    for (size_t i = 1; i < lines.size(); ++i) {
        std::string line = PSFParserStringUtils::trim(lines[i]);
        // Skip lines that look like section headers
        if (line.find('!') != std::string::npos) {
            continue;
        }

        std::istringstream iss_line(line);
        int idx;
        while (iss_line >> idx) {
            // Skip zero values as they are placeholders in PSF format
            if (idx == 0) {
                continue;
            }
            // Convert to 0-based indexing and check validity
            idx--;  // Convert to 0-based indexing
            if (idx < 0 || idx >= topology.get_num_atoms()) {
                std::cerr << "Invalid atom index " << (idx + 1) << " in bond" << std::endl;
                return false;
            }
            bond_indices.push_back(idx);
        }
    }

    // Add bonds to topology
    for (size_t i = 0; i < bond_indices.size(); i += 2) {
        topology.add_bond(bond_indices[i], bond_indices[i + 1]);
    }

    return true;
}

bool PSFParserSectionsBasic::parse_angles_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    // Parse number of angles from the first line
    std::istringstream iss(PSFParserStringUtils::trim(lines[0]));
    int num_angles = 0;
    if (!(iss >> num_angles)) {
        std::cerr << "Failed to parse number of angles" << std::endl;
        return false;
    }

    if (num_angles <= 0) {
        // Angles are optional, so just return true if there are none
        return true;
    }

    // Read angle indices
    std::vector<int> angle_indices;
    angle_indices.reserve(num_angles * 3);

    for (size_t i = 1; i < lines.size(); ++i) {
        std::string line = PSFParserStringUtils::trim(lines[i]);
        std::istringstream iss(line);
        int idx;
        while (iss >> idx) {
            // Skip zero values as they are placeholders in PSF format
            if (idx == 0) {
                continue;
            }
            // Convert to 0-based indexing and check validity
            idx--;  // Convert to 0-based indexing
            if (idx < 0 || idx >= topology.get_num_atoms()) {
                std::cerr << "Invalid atom index " << (idx + 1) << " in angle" << std::endl;
                return false;
            }
            angle_indices.push_back(idx);
        }
    }

    // Add angles to topology
    for (size_t i = 0; i < angle_indices.size(); i += 3) {
        topology.add_angle(
            angle_indices[i],
            angle_indices[i + 1],
            angle_indices[i + 2]
        );
    }

    return true;
}

bool PSFParserSectionsBasic::parse_dihedrals_from_lines(const std::vector<std::string>& dihedral_lines,
                                                       model::Topology& topology, 
                                                       const std::string& /* section_name */) {
    if (dihedral_lines.empty()) {
        std::cerr << "No lines provided for dihedrals section." << std::endl;
        return false;
    }

    // The first line should contain the number of dihedrals
    std::istringstream iss(PSFParserStringUtils::trim(dihedral_lines[0]));
    int num_dihedrals = 0;
    if (!(iss >> num_dihedrals)) {
        std::cerr << "Failed to parse dihedral count from line: " << dihedral_lines[0] << std::endl;
        return false;
    }

    // If there's no dihedral to parse, nothing to do
    if (num_dihedrals <= 0) {
        return true;
    }

    // First collect all dihedral indices
    std::vector<int> all_indices;
    all_indices.reserve(num_dihedrals * 4);

    // Begin from line 1 because line 0 is the count
    for (size_t line_idx = 1; line_idx < dihedral_lines.size(); ++line_idx) {
        std::istringstream iss_line(PSFParserStringUtils::trim(dihedral_lines[line_idx]));
        int atom_idx;
        while (iss_line >> atom_idx) {
            if (atom_idx == 0) continue;  // Skip fillers
            atom_idx--;  // Convert to 0-based
            if (atom_idx < 0 || atom_idx >= topology.get_num_atoms()) {
                std::cerr << "Invalid atom index in dihedral: " << (atom_idx + 1) << std::endl;
                return false;
            }
            all_indices.push_back(atom_idx);
        }
    }

    // Check if we got the expected number of indices
    const size_t expected_indices = static_cast<size_t>(num_dihedrals) * 4;
    if (all_indices.size() != expected_indices) {
        std::cerr << "Wrong number of dihedral indices: expected " << expected_indices
                  << ", got " << all_indices.size() << std::endl;
        return false;
    }

    // Add all dihedrals to topology
    for (size_t i = 0; i < all_indices.size(); i += 4) {
        topology.add_dihedral(all_indices[i], all_indices[i+1], 
                            all_indices[i+2], all_indices[i+3]);
    }

    return true;
}

bool PSFParserSectionsBasic::parse_dihedrals_from_lines(const std::vector<std::string>& dihedral_lines,
                                                       model::Topology& topology) {
    return parse_dihedrals_from_lines(dihedral_lines, topology, "");
}

} // namespace io
} // namespace pygcmc