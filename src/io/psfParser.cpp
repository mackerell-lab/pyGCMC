#include "psfParser.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <unordered_map>

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

    // Define section order for consistent parsing
    const std::vector<std::string> sections_order = {
        "NATOM", "NBOND", "NTHETA", "NPHI", "NIMPHI", "CMAP", "NGRP", "NDON", "NACC"
    };

    // Find start of each section
    std::unordered_map<std::string, size_t> section_index;
    for (size_t i = 0; i < lines.size(); i++) {
        std::string trimmed = trim(lines[i]);
        if (trimmed.find("!NATOM") != std::string::npos) {
            section_index["NATOM"] = i;
        } else if (trimmed.find("!NBOND") != std::string::npos) {
            section_index["NBOND"] = i;
        } else if (trimmed.find("!NTHETA") != std::string::npos) {
            section_index["NTHETA"] = i;
        } else if (trimmed.find("!NPHI") != std::string::npos) {
            section_index["NPHI"] = i;
        } else if (trimmed.find("!NIMPHI") != std::string::npos) {
            section_index["NIMPHI"] = i;
        } else if (trimmed.find("!NCRTERM") != std::string::npos ||
                   trimmed.find("!NCMAP") != std::string::npos) {
            section_index["CMAP"] = i;
        } else if (trimmed.find("!NGRP") != std::string::npos) {
            section_index["NGRP"] = i;
        } else if (trimmed.find("!NDON") != std::string::npos) {
            section_index["NDON"] = i;
        } else if (trimmed.find("!NACC") != std::string::npos) {
            section_index["NACC"] = i;
        }
    }

    // Extract each section's lines into dedicated buffers
    std::unordered_map<std::string, std::vector<std::string>> section_data;
    for (size_t s = 0; s < sections_order.size(); s++) {
        const std::string& key = sections_order[s];
        if (!section_index.count(key)) continue;

        size_t start_line = section_index[key];
        
        // Find next section's start or end of file
        size_t end_line = lines.size();
        for (size_t t = s + 1; t < sections_order.size(); t++) {
            const std::string& next_key = sections_order[t];
            if (section_index.count(next_key)) {
                size_t next_start = section_index[next_key];
                if (next_start > start_line && next_start < end_line) {
                    end_line = next_start;
                }
            }
        }

        // Validate section size
        if (start_line >= end_line || start_line >= lines.size()) {
            std::cerr << "Invalid section bounds for " << key << std::endl;
            return false;
        }

        // Copy lines for this section
        section_data[key].reserve(end_line - start_line);
        for (size_t l = start_line; l < end_line; l++) {
            section_data[key].push_back(lines[l]);
        }

        // Validate we have at least one line (for count)
        if (section_data[key].empty()) {
            std::cerr << "Empty section: " << key << std::endl;
            return false;
        }
    }

    // Parse sections in canonical order
    size_t current_line = 0;  // Each section starts at line 0

    // NATOM must be first
    if (!section_data.count("NATOM")) {
        std::cerr << "No NATOM section found" << std::endl;
        return false;
    }

    // Reserve space for atoms to avoid reallocation
    std::istringstream iss(trim(section_data["NATOM"][0]));
    int num_atoms = 0;
    if (iss >> num_atoms && num_atoms > 0) {
        topology.reserve_atoms(num_atoms);
    }

    if (!parse_atoms_from_lines(section_data["NATOM"], current_line, topology)) {
        return false;
    }

    // Parse optional sections in order
    if (section_data.count("NBOND")) {
        current_line = 0;  // Reset for new section
        if (!parse_bonds_from_lines(section_data["NBOND"], current_line, topology)) return false;
    }

    if (section_data.count("NTHETA")) {
        current_line = 0;  // Reset for new section
        if (!parse_angles_from_lines(section_data["NTHETA"], current_line, topology)) return false;
    }

    if (section_data.count("NPHI")) {
        current_line = 0;  // Reset for new section
        if (!parse_dihedrals_from_lines(section_data["NPHI"], current_line, topology)) return false;
    }

    if (section_data.count("NIMPHI")) {
        current_line = 0;  // Reset for new section
        if (!parse_impropers_from_lines(section_data["NIMPHI"], current_line, topology)) return false;
    }

    if (section_data.count("CMAP")) {
        current_line = 0;  // Reset for new section
        if (!parse_cmap_from_lines(section_data["CMAP"], current_line, topology)) return false;
    }

    if (section_data.count("NGRP")) {
        current_line = 0;  // Reset for new section
        if (!parse_groups_from_lines(section_data["NGRP"], current_line, topology)) return false;
    }

    if (section_data.count("NDON")) {
        current_line = 0;  // Reset for new section
        if (!parse_donors_from_lines(section_data["NDON"], current_line, topology)) return false;
    }

    if (section_data.count("NACC")) {
        current_line = 0;  // Reset for new section
        if (!parse_acceptors_from_lines(section_data["NACC"], current_line, topology)) return false;
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

bool PSFParser::parse_atoms_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology) {
    // Parse number of atoms from the current line
    std::istringstream iss(trim(lines[current_line]));
    int num_atoms = 0;
    if (!(iss >> num_atoms)) {
        std::cerr << "Failed to parse number of atoms" << std::endl;
        return false;
    }

    if (num_atoms <= 0) {
        std::cerr << "Invalid number of atoms: " << num_atoms << std::endl;
        return false;
    }

    // Read exactly num_atoms lines
    for (int i = 0; i < num_atoms; ++i) {
        current_line++;
        if (current_line >= lines.size()) {
            std::cerr << "Failed to read atom line " << i + 1 << std::endl;
            return false;
        }

        std::istringstream iss(lines[current_line]);
        int atomIndex;
        std::string segment_name;
        int residue_number;
        std::string residue_name;
        std::string atom_name;
        std::string atom_type;
        double charge;
        double mass;
        int unusedField;

        bool parsedOK = false;
        if (iss >> atomIndex >> segment_name >> residue_number
                >> residue_name >> atom_name >> atom_type
                >> charge >> mass >> unusedField) {
            parsedOK = true;
        } else {
            iss.clear();
            iss.seekg(0);
            if (iss >> atomIndex >> segment_name >> residue_number
                    >> residue_name >> atom_name >> atom_type
                    >> charge >> mass) {
                parsedOK = true;
            }
        }

        if (!parsedOK) {
            std::cerr << "Failed to parse atom line " << (i + 1)
                      << ": " << lines[current_line] << std::endl;
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
    }

    return true;
}

bool PSFParser::parse_bonds_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology) {
    // Parse number of bonds from the current line
    std::istringstream iss(trim(lines[current_line]));
    int num_bonds = 0;
    if (!(iss >> num_bonds)) {
        std::cerr << "Failed to parse number of bonds" << std::endl;
        return false;
    }

    if (num_bonds <= 0) {
        // Bonds are optional, so just return true if there are none
        return true;
    }

    // Read bond indices
    std::vector<int> bond_indices;
    bond_indices.reserve(num_bonds * 2);

    while (bond_indices.size() < static_cast<size_t>(num_bonds) * 2) {
        current_line++;
        if (current_line >= lines.size()) {
            std::cerr << "Unexpected end of file while reading bonds" << std::endl;
            return false;
        }

        std::istringstream iss(lines[current_line]);
        int idx;
        while (iss >> idx) {
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

bool PSFParser::parse_angles_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology) {
    // Parse number of angles from the current line
    std::istringstream iss(trim(lines[current_line]));
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

    while (angle_indices.size() < static_cast<size_t>(num_angles) * 3) {
        current_line++;
        if (current_line >= lines.size()) {
            std::cerr << "Unexpected end of file while reading angles" << std::endl;
            return false;
        }

        std::istringstream iss(lines[current_line]);
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

bool PSFParser::parse_dihedrals_from_lines(const std::vector<std::string>& lines,
                                         size_t& current_line,
                                         model::Topology& topology) {
    if (current_line >= lines.size()) {
        std::cerr << "Invalid line number in parse_dihedrals_from_lines" << std::endl;
        return false;
    }

    // Parse number of dihedrals from the current line
    std::istringstream iss(trim(lines[current_line]));
    size_t num_dihedrals = 0;
    if (!(iss >> num_dihedrals)) {
        std::cerr << "Failed to parse number of dihedrals" << std::endl;
        return false;
    }

    if (num_dihedrals <= 0) {
        return true;  // Dihedrals are optional
    }

    size_t dihedrals_parsed = 0;
    std::vector<int> buffer;
    buffer.reserve(4);  // Each dihedral has 4 atoms

    // Calculate how many lines we need (assuming 8 indices per line = 2 dihedrals)
    size_t min_lines_needed = (num_dihedrals * 4 + 7) / 8;  // Round up division
    if (current_line + 1 + min_lines_needed > lines.size()) {
        std::cerr << "Not enough lines for dihedrals: need " << min_lines_needed 
                  << ", have " << (lines.size() - current_line - 1) << std::endl;
        return false;
    }

    for (size_t i = current_line + 1; i < lines.size() && dihedrals_parsed < num_dihedrals; i++) {
        std::istringstream iss_line(trim(lines[i]));
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) continue; // Skip placeholder zeros
            idx--;  // Convert to 0-based indexing
            if (idx < 0 || idx >= topology.get_num_atoms()) {
                std::cerr << "Invalid atom index " << (idx + 1) << " in dihedral" << std::endl;
                return false;
            }
            buffer.push_back(idx);
            if (buffer.size() == 4) {
                topology.add_dihedral(buffer[0], buffer[1], buffer[2], buffer[3], 1, 0.0, 0.0, false);
                buffer.clear();
                dihedrals_parsed++;
                if (dihedrals_parsed == num_dihedrals) {
                    break;
                }
            }
        }
    }

    if (dihedrals_parsed != num_dihedrals) {
        std::cerr << "Failed to parse all dihedrals: expected " << num_dihedrals 
                  << ", got " << dihedrals_parsed << std::endl;
        return false;
    }

    return true;
}

bool PSFParser::parse_impropers_from_lines(const std::vector<std::string>& lines,
                                         size_t& current_line,
                                         model::Topology& topology) {
    if (current_line >= lines.size()) {
        std::cerr << "Invalid line number in parse_impropers_from_lines" << std::endl;
        return false;
    }

    std::istringstream iss(trim(lines[current_line]));
    size_t num_impropers = 0;
    if (!(iss >> num_impropers)) {
        std::cerr << "Failed to parse number of impropers" << std::endl;
        return false;
    }

    if (num_impropers <= 0) {
        return true;  // Impropers are optional
    }

    size_t impropers_parsed = 0;
    std::vector<int> buffer;
    buffer.reserve(4);  // Each improper has 4 atoms

    // Calculate how many lines we need (assuming 8 indices per line = 2 impropers)
    size_t min_lines_needed = (num_impropers * 4 + 7) / 8;  // Round up division
    if (current_line + 1 + min_lines_needed > lines.size()) {
        std::cerr << "Not enough lines for impropers: need " << min_lines_needed 
                  << ", have " << (lines.size() - current_line - 1) << std::endl;
        return false;
    }

    for (size_t i = current_line + 1; i < lines.size() && impropers_parsed < num_impropers; i++) {
        std::istringstream iss_line(trim(lines[i]));
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) continue; // Skip placeholder zeros
            idx--;  // Convert to 0-based indexing
            if (idx < 0 || idx >= topology.get_num_atoms()) {
                std::cerr << "Invalid atom index " << (idx + 1) << " in improper" << std::endl;
                return false;
            }
            buffer.push_back(idx);
            if (buffer.size() == 4) {
                topology.add_improper(buffer[0], buffer[1], buffer[2], buffer[3]);
                buffer.clear();
                impropers_parsed++;
                if (impropers_parsed == num_impropers) {
                    break;
                }
            }
        }
    }

    if (impropers_parsed != num_impropers) {
        std::cerr << "Failed to parse all impropers: expected " << num_impropers 
                  << ", got " << impropers_parsed << std::endl;
        return false;
    }

    return true;
}

bool PSFParser::parse_donors_from_lines(const std::vector<std::string>& lines,
                                      size_t& current_line,
                                      model::Topology& topology) {
    std::istringstream iss(trim(lines[current_line]));
    size_t num_donors = 0;
    if (!(iss >> num_donors)) {
        std::cerr << "Failed to parse number of donors" << std::endl;
        return false;
    }

    size_t donors_parsed = 0;
    for (size_t i = current_line + 1; i < lines.size() && donors_parsed < num_donors; i++) {
        std::istringstream iss_line(trim(lines[i]));
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

bool PSFParser::parse_acceptors_from_lines(const std::vector<std::string>& lines,
                                         size_t& current_line,
                                         model::Topology& topology) {
    std::istringstream iss(trim(lines[current_line]));
    size_t num_acceptors = 0;
    if (!(iss >> num_acceptors)) {
        std::cerr << "Failed to parse number of acceptors" << std::endl;
        return false;
    }

    size_t acceptors_parsed = 0;
    for (size_t i = current_line + 1; i < lines.size() && acceptors_parsed < num_acceptors; i++) {
        std::istringstream iss_line(trim(lines[i]));
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

bool PSFParser::parse_cmap_from_lines(const std::vector<std::string>& lines,
                                    size_t& current_line,
                                    model::Topology& topology) {
    std::istringstream iss(trim(lines[current_line]));
    size_t num_cmaps = 0;
    if (!(iss >> num_cmaps)) {
        std::cerr << "Failed to parse number of CMAP terms" << std::endl;
        return false;
    }

    size_t cmaps_parsed = 0;
    std::vector<int> buffer;
    for (size_t i = current_line + 1; i < lines.size() && cmaps_parsed < num_cmaps; i++) {
        std::istringstream iss_line(trim(lines[i]));
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

bool PSFParser::parse_groups_from_lines(const std::vector<std::string>& lines,
                                      size_t& current_line,
                                      model::Topology& topology) {
    std::istringstream iss(trim(lines[current_line]));
    size_t num_groups = 0;
    if (!(iss >> num_groups)) {
        std::cerr << "Failed to parse number of groups" << std::endl;
        return false;
    }

    size_t groups_parsed = 0;
    std::vector<int> current_group;
    int current_group_id = 1;  // Start with group ID 1

    for (size_t i = current_line + 1; i < lines.size() && groups_parsed < num_groups; i++) {
        std::istringstream iss_line(trim(lines[i]));
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
