#include "psfParser.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <set>

namespace pygcmc {
namespace io {

// Helper function to read file lines
bool readFileToLines(const std::string& filename, std::vector<std::string>& lines) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    return true;
}

std::string PSFParser::trim(const std::string& str) {
    const auto start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    const auto end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

bool PSFParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    // Read all lines at once
    std::vector<std::string> lines;
    if (!readFileToLines(filename, lines)) {
        return false;
    }

    // First, collect all sections
    std::unordered_map<std::string, std::vector<std::string>> sections;
    std::string current_section;
    std::vector<std::string> current_lines;

    for (const auto& line : lines) {
        std::string trimmed = trim(line);
        if (trimmed.empty()) continue;

        // Check if this is a section header
        if (trimmed.find('!') != std::string::npos) {
            // Save previous section if any
            if (!current_section.empty() && !current_lines.empty()) {
                sections[current_section] = current_lines;
            }

            // Determine new section type
            if (trimmed.find("!NATOM") != std::string::npos) {
                current_section = "NATOM";
            } else if (trimmed.find("!NBOND") != std::string::npos) {
                current_section = "NBOND";
            } else if (trimmed.find("!NTHETA") != std::string::npos) {
                current_section = "NTHETA";
            } else if (trimmed.find("!NPHI") != std::string::npos) {
                current_section = "NPHI";
            } else if (trimmed.find("!NIMPHI") != std::string::npos) {
                current_section = "NIMPHI";
            } else if (trimmed.find("!NCRTERM") != std::string::npos ||
                      trimmed.find("!NCMAP") != std::string::npos) {
                current_section = "CMAP";
            } else if (trimmed.find("!NGRP") != std::string::npos) {
                current_section = "NGRP";
            } else if (trimmed.find("!NDON") != std::string::npos) {
                current_section = "NDON";
            } else if (trimmed.find("!NACC") != std::string::npos) {
                current_section = "NACC";
            } else {
                current_section = "";  // Unknown section
            }
            current_lines.clear();
            current_lines.push_back(trimmed);  // Include the header line
        } else if (!current_section.empty()) {
            // Add line to current section
            current_lines.push_back(trimmed);
        }
    }

    // Add the last section if any
    if (!current_section.empty() && !current_lines.empty()) {
        sections[current_section] = current_lines;
    }

    // Now process sections in the required order
    // NATOM must be first
    if (!sections.count("NATOM")) {
        std::cerr << "Missing required NATOM section" << std::endl;
        return false;
    }

    if (!parse_atoms_from_lines(sections["NATOM"], topology)) {
        std::cerr << "Failed to parse NATOM section" << std::endl;
        return false;
    }

    // Process optional sections in order
    const std::vector<std::string> optional_sections = {
        "NBOND", "NTHETA", "NPHI", "NIMPHI", "CMAP", "NGRP", "NDON", "NACC"
    };

    for (const auto& section : optional_sections) {
        if (sections.count(section)) {
            bool success = false;

            if (section == "NBOND") {
                success = parse_bonds_from_lines(sections[section], topology);
            } else if (section == "NTHETA") {
                success = parse_angles_from_lines(sections[section], topology);
            } else if (section == "NPHI") {
                success = parse_dihedrals_from_lines(sections[section], topology);
            } else if (section == "NIMPHI") {
                success = parse_impropers_from_lines(sections[section], topology);
            } else if (section == "CMAP") {
                success = parse_cmap_from_lines(sections[section], topology);
            } else if (section == "NGRP") {
                success = parse_groups_from_lines(sections[section], topology);
            } else if (section == "NDON") {
                success = parse_donors_from_lines(sections[section], topology);
            } else if (section == "NACC") {
                success = parse_acceptors_from_lines(sections[section], topology);
            }

            if (!success) {
                std::cerr << "Failed to parse " << section << " section" << std::endl;
                return false;
            }
        }
    }

    return true;
}

bool PSFParser::parse_title_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology) {
    // Parse number of titles from the current line
    std::istringstream iss(trim(lines[current_line]));
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

bool PSFParser::parse_atoms_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    // Parse number of atoms from the first line (header)
    std::istringstream iss(trim(lines[0]));  // Always use first line as header
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
        std::string line = trim(lines[i]);
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

bool PSFParser::parse_bonds_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    // Parse number of bonds from the first line
    std::istringstream iss(trim(lines[0]));
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
        std::string line = trim(lines[i]);
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

bool PSFParser::parse_angles_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    // Parse number of angles from the first line
    std::istringstream iss(trim(lines[0]));
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
        std::string line = trim(lines[i]);
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

bool PSFParser::parse_dihedrals_from_lines(const std::vector<std::string>& dihedral_lines,
                                           model::Topology& topology)
{
    if (dihedral_lines.empty()) {
        std::cerr << "No lines provided for dihedrals section." << std::endl;
        return false;
    }

    // The first line should contain the number of dihedrals
    std::istringstream iss(trim(dihedral_lines[0]));
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
        std::istringstream iss_line(trim(dihedral_lines[line_idx]));
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

bool PSFParser::parse_impropers_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    int num_impropers = 0;
    if (!(iss >> num_impropers)) {
        std::cerr << "Failed to parse number of impropers" << std::endl;
        return false;
    }

    if (num_impropers <= 0) {
        // no impropers
        return true;
    }

    int impropers_added = 0;
    std::vector<int> tmp_indices;

    for (size_t i = 1; i < lines.size(); ++i) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) continue;
            idx--; // 0-based
            if (idx < 0 || idx >= topology.get_num_atoms()) {
                std::cerr << "Invalid atom index in impropers" << std::endl;
                return false;
            }
            tmp_indices.push_back(idx);
            if (tmp_indices.size() == 4) {
                topology.add_improper(tmp_indices[0],
                                    tmp_indices[1],
                                    tmp_indices[2],
                                    tmp_indices[3]);
                tmp_indices.clear();
                impropers_added++;
            }
        }
    }

    return impropers_added == num_impropers;
}

bool PSFParser::parse_donors_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_donors = 0;
    if (!(iss >> num_donors)) {
        std::cerr << "Failed to parse number of donors" << std::endl;
        return false;
    }

    size_t donors_parsed = 0;
    for (size_t i = 1; i < lines.size() && donors_parsed < num_donors; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int donor_idx, hydrogen_idx;
        while (iss_line >> donor_idx >> hydrogen_idx) {
            if (donor_idx == 0 || hydrogen_idx == 0) continue;
            donor_idx--; // Convert to 0-based indexing
            hydrogen_idx--;
            topology.add_donor(donor_idx, hydrogen_idx);
            donors_parsed++;
            if (donors_parsed == num_donors) {
                break;
            }
        }
    }
    return donors_parsed == num_donors;
}

bool PSFParser::parse_acceptors_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_acceptors = 0;
    if (!(iss >> num_acceptors)) {
        std::cerr << "Failed to parse number of acceptors" << std::endl;
        return false;
    }

    size_t acceptors_parsed = 0;
    for (size_t i = 1; i < lines.size() && acceptors_parsed < num_acceptors; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int acceptor_idx;
        while (iss_line >> acceptor_idx) {
            if (acceptor_idx == 0) continue;
            acceptor_idx--; // Convert to 0-based indexing
            topology.add_acceptor(acceptor_idx);
            acceptors_parsed++;
            if (acceptors_parsed == num_acceptors) {
                break;
            }
        }
    }
    return acceptors_parsed == num_acceptors;
}

bool PSFParser::parse_cmap_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_cmaps = 0;
    if (!(iss >> num_cmaps)) {
        std::cerr << "Failed to parse number of CMAP terms" << std::endl;
        return false;
    }

    size_t cmaps_parsed = 0;
    std::vector<int> buffer;
    for (size_t i = 1; i < lines.size() && cmaps_parsed < num_cmaps; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) continue; // Skip placeholder zeros
            buffer.push_back(idx - 1); // Convert to 0-based indexing
            if (buffer.size() == 8) {
                std::array<int, 8> cmap_array;
                std::copy_n(buffer.begin(), 8, cmap_array.begin());
                topology.add_cmap(cmap_array);
                buffer.clear();
                cmaps_parsed++;
                if (cmaps_parsed == num_cmaps) {
                    break;
                }
            }
        }
    }
    return cmaps_parsed == num_cmaps;
}

bool PSFParser::parse_groups_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_groups = 0;
    if (!(iss >> num_groups)) {
        std::cerr << "Failed to parse number of groups" << std::endl;
        return false;
    }

    size_t groups_parsed = 0;
    std::vector<int> current_group;
    int current_group_id = 1;  // Start with group ID 1

    for (size_t i = 1; i < lines.size() && groups_parsed < num_groups; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) {
                // Zero marks end of current group
                if (!current_group.empty()) {
                    topology.add_group(current_group_id++, current_group, "");  // Empty type string
                    current_group.clear();
                    groups_parsed++;
                    if (groups_parsed == num_groups) {
                        break;
                    }
                }
                continue;
            }
            current_group.push_back(idx - 1); // Convert to 0-based indexing
        }
    }
    // Add last group if not empty
    if (!current_group.empty()) {
        topology.add_group(current_group_id, current_group, "");
        groups_parsed++;
    }
    return groups_parsed == num_groups;
}

} // namespace io
} // namespace pygcmc
