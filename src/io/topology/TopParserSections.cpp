// src/io/topology/TopParserSections.cpp

#include "TopParserSections.hpp"
#include <algorithm>
#include <array>

namespace pygcmc {
namespace io {

bool TopParserSections::parse_moleculetype_section(const std::vector<LineInfo>& lines, [[maybe_unused]] model::Topology& topology,
                                                  std::string& current_molecule_type, int& current_molecule_nrexcl) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 2) continue;

        try {
            current_molecule_type = tokens[0];
            current_molecule_nrexcl = std::stoi(tokens[1]);
            return true;
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Error parsing moleculetype line: ", line, "\n");
            return false;
        }
    }
    return true;
}

bool TopParserSections::parse_atoms_section(const std::vector<LineInfo>& lines, model::Topology& topology,
                                           const std::string& current_molecule_type) {
    // Get current residue number offset from segment name
    int segment_index = 0;
    size_t underscore_pos = current_molecule_type.find_last_of('_');
    if (underscore_pos != std::string::npos) {
        try {
            segment_index = std::stoi(current_molecule_type.substr(underscore_pos + 1));
        } catch (...) {
            segment_index = 0;
        }
    }

    // Find the maximum residue number in this molecule definition
    int max_resnum = 0;
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 7) continue;
        try {
            int resnum = std::stoi(tokens[2]);
            max_resnum = std::max(max_resnum, resnum);
        } catch (...) {
            continue;
        }
    }

    // Now parse atoms with adjusted residue numbers
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 7) continue;  // Need at least 7 columns (tip3p.itp format)

        try {
            // Parse atom data
            std::string atom_type = tokens[1];
            int residue_number = std::stoi(tokens[2]);
            std::string residue_name = tokens[3];
            std::string atom_name = tokens[4];
            double charge = std::stod(tokens[6]);
            // Use explicit mass if present; otherwise infer default
            double mass = (tokens.size() >= 8) ? std::stod(tokens[7])
                                              : TopParserUtilities::default_mass_for_atom_type(atom_type);

            // Adjust residue number based on segment index
            int adjusted_resnum = residue_number + (segment_index * max_resnum);

            // Add atom to topology using current_molecule_type as segment
            topology.add_atom(
                atom_name,
                atom_type,
                charge,
                mass,
                residue_name,
                adjusted_resnum,
                current_molecule_type
            );
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Error parsing atom line: ", line, "\n");
            return false;
        }
    }
    return true;
}

bool TopParserSections::parse_bonds_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 2) continue;  // Need at least ai, aj

        try {
            // Always parse the first two atoms (required)
            int atom1 = std::stoi(tokens[0]) - 1 + atom_offset;  // Convert to 0-based indexing and add offset
            int atom2 = std::stoi(tokens[1]) - 1 + atom_offset;

            // Parse function type if present (default to 1)
            int func_type = 1;
            if (tokens.size() >= 3) {
                func_type = std::stoi(tokens[2]);
            }

            // Parse first set of parameters if present
            double length = 0.0;
            double force_const = 0.0;
            if (tokens.size() >= 5) {
                length = std::stod(tokens[3]);
                force_const = std::stod(tokens[4]);
            }

            // Note: We ignore any additional parameter sets (tokens[5] onwards)
            // as they are typically alternative parameters for different force field variants

            // Add bond to topology
            topology.add_bond(atom1, atom2, length, force_const, func_type);
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Warning: Error parsing bond line: ", line, " at ",
                     line_info.source_file, ":", line_info.line_number,
                     " - ", e.what(), "\n");
            continue;  // Continue with next line instead of failing
        }
    }
    return true;  // Return true if we processed all lines (even with some warnings)
}

bool TopParserSections::parse_angles_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 3) continue;  // Need at least ai, aj, ak

        try {
            // Always parse the three atoms (required)
            int atom1 = std::stoi(tokens[0]) - 1 + atom_offset;
            int atom2 = std::stoi(tokens[1]) - 1 + atom_offset;
            int atom3 = std::stoi(tokens[2]) - 1 + atom_offset;

            // Parse function type if present (default to 1)
            int func_type = 1;
            if (tokens.size() >= 4) {
                func_type = std::stoi(tokens[3]);
            }

            // Parse first set of parameters if present
            double angle = 0.0;
            double force_const = 0.0;
            if (tokens.size() >= 6) {
                angle = std::stod(tokens[4]);
                force_const = std::stod(tokens[5]);
            }

            // Note: We ignore any additional parameter sets (tokens[6] onwards)
            // as they are typically alternative parameters for different force field variants

            // Add angle to topology
            topology.add_angle(atom1, atom2, atom3, angle, force_const, func_type);
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Warning: Error parsing angle line: ", line, " at ",
                     line_info.source_file, ":", line_info.line_number,
                     " - ", e.what(), "\n");
            continue;  // Continue with next line instead of failing
        }
    }
    return true;
}

bool TopParserSections::parse_dihedrals_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    TopParserUtilities::debug_print("Parsing ", lines.size(), " dihedral lines\n");
    int proper_count = 0;
    int improper_count = 0;

    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 4) continue;

        try {
            int atom1 = std::stoi(tokens[0]) - 1 + atom_offset;
            int atom2 = std::stoi(tokens[1]) - 1 + atom_offset;
            int atom3 = std::stoi(tokens[2]) - 1 + atom_offset;
            int atom4 = std::stoi(tokens[3]) - 1 + atom_offset;

            // Default to function type 1 (proper dihedral) if not specified
            int funcType = (tokens.size() >= 5) ? std::stoi(tokens[4]) : 1;

            if (funcType == 2 || funcType == 4) {
                topology.add_improper(atom1, atom2, atom3, atom4);
                improper_count++;
            } else {
                topology.add_dihedral(atom1, atom2, atom3, atom4);
                proper_count++;
                if (funcType == 9 && tokens.size() >= 7) {
                    int multiplicity = std::stoi(tokens[6]);
                    for (int i = 1; i < multiplicity; i++) {
                        topology.add_dihedral(atom1, atom2, atom3, atom4);
                        proper_count++;
                    }
                }
            }
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Error parsing dihedral line: ", line, " - ", e.what(), "\n");
            return false;
        }
    }

    TopParserUtilities::debug_print("Added ", proper_count, " proper dihedrals and ",
              improper_count, " improper dihedrals\n");
    return true;
}

bool TopParserSections::parse_impropers_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 4) continue;  // Need at least ai, aj, ak, al

        try {
            int atom1 = std::stoi(tokens[0]) - 1 + atom_offset;
            int atom2 = std::stoi(tokens[1]) - 1 + atom_offset;
            int atom3 = std::stoi(tokens[2]) - 1 + atom_offset;
            int atom4 = std::stoi(tokens[3]) - 1 + atom_offset;

            // Add improper to topology
            topology.add_improper(atom1, atom2, atom3, atom4);
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Error parsing improper line: ", line, "\n");
            return false;
        }
    }
    return true;
}

bool TopParserSections::parse_molecules_section(const std::vector<LineInfo>& lines, [[maybe_unused]] model::Topology& topology,
                                               std::vector<std::pair<std::string, int>>& molecule_order) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 2) continue;

        try {
            std::string mol_type = tokens[0];
            int count = std::stoi(tokens[1]);

            // Store molecule type and count
            molecule_order.push_back(std::make_pair(mol_type, count));
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Error parsing molecules line: ", line, "\n");
            return false;
        }
    }
    return true;
}

bool TopParserSections::parse_cmaps_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line));
        if (tokens.size() < 6) {  // Need 5 atoms + function type
            TopParserUtilities::debug_print("Warning: Skipping CMAP line with insufficient tokens: ", line, "\n");
            continue;
        }

        try {
            // GROMACS format: ai aj ak al am funct
            // These atoms define the two consecutive phi-psi dihedrals
            std::array<int, 5> cmap_atoms;
            for (int i = 0; i < 5; ++i) {
                cmap_atoms[i] = std::stoi(tokens[i]) - 1 + atom_offset;  // Convert to 0-based indexing
            }
            int function_type = std::stoi(tokens[5]);

            // Add CMAP to topology using the GROMACS format overload
            topology.add_cmap(cmap_atoms, function_type);

            TopParserUtilities::debug_print("Added CMAP between atoms: ", cmap_atoms[0], " ", cmap_atoms[1], " ", cmap_atoms[2], " ", cmap_atoms[3], " ", cmap_atoms[4], " (function type ", function_type, ")\n");
        } catch (const std::exception& e) {
            TopParserUtilities::debug_print("Error parsing CMAP line: ", line, " at ",
                     line_info.source_file, ":", line_info.line_number,
                     " - ", e.what(), "\n");
            return false;
        }
    }
    return true;
}

} // namespace io
} // namespace pygcmc
