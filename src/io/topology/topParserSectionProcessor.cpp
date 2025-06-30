// src/io/topology/topParserSectionProcessor.cpp

#include "topParserSectionProcessor.hpp"
#include "topParserUtilities.hpp"
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace io {

SectionProcessorResult TopParserSectionProcessor::process_sections(const std::vector<LineInfo>& all_lines) {
    SectionProcessorResult result;
    
    // Parse sections using state machine
    std::string current_mol_type;
    bool inside_molecule = false;
    std::string current_section;

    for (const auto& line_info : all_lines) {
        std::string line_no_comment = TopParserUtilities::remove_comment(line_info.content);
        std::string trimmed = TopParserUtilities::trim(line_no_comment);
        if (trimmed.empty()) {
            continue;
        }

        // Check if line starts with '['
        if (trimmed.size() > 1 && trimmed[0] == '[' && trimmed[trimmed.size() - 1] == ']') {
            // Remove the leading '[' and trailing ']' and trim whitespace
            std::string section_name = trimmed.substr(1, trimmed.size() - 2);
            section_name = TopParserUtilities::trim(section_name);

            // Convert to lowercase
            std::transform(section_name.begin(), section_name.end(), 
                         section_name.begin(), ::tolower);

            if (section_name == "moleculetype") {
                inside_molecule = true;
                current_mol_type.clear();
                current_section = "moleculetype";
                result.found_valid_section = true;
                TopParserUtilities::debug_print("\n=== Found [ moleculetype ] section at ", line_info.source_file, 
                         ":" , line_info.line_number, " ===\n");
            }
            else if (section_name == "molecules") {
                inside_molecule = false;
                current_mol_type.clear();
                current_section = "molecules";
                result.found_valid_section = true;
                TopParserUtilities::debug_print("\n=== Found [ molecules ] section at ", line_info.source_file, 
                         ":" , line_info.line_number, " ===\n");
            }
            else if (inside_molecule) {
                current_section = section_name;
                result.found_valid_section = true;
                TopParserUtilities::debug_print("Found [ ", section_name, " ] section for molecule ", 
                         current_mol_type, "\n");
            }
            else {
                current_section.clear();
            }
            continue;
        }

        // Process content based on current state
        if (inside_molecule) {
            if (current_section == "moleculetype") {
                process_moleculetype_section(trimmed, line_info,
                                            current_mol_type, 
                                            result.molecule_types_order,
                                            result.molecule_atoms_temp,
                                            result.molecule_bonds_temp,
                                            result.molecule_angles_temp,
                                            result.molecule_dihedrals_temp,
                                            result.molecule_impropers_temp);
            }
            else if (!current_mol_type.empty()) {
                if (current_section == "atoms") {
                    result.molecule_atoms_temp[current_mol_type].push_back(line_info);
                    auto tokens = TopParserUtilities::split(trimmed);
                    if (tokens.size() >= 4) {
                        TopParserUtilities::debug_print("Added atom to ", current_mol_type, ": ", 
                                 tokens[0], " ", tokens[1], " (residue ", 
                                 tokens[3], ")\n");
                    }
                }
                else if (current_section == "settles") {
                    TopParserUtilities::debug_print("Found settles for ", current_mol_type, ": ", 
                             trimmed, "\n");
                }
                else if (current_section == "bonds") {
                    result.molecule_bonds_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "angles") {
                    result.molecule_angles_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "dihedrals") {
                    result.molecule_dihedrals_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "impropers") {
                    result.molecule_impropers_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "cmap") {
                    result.molecule_cmaps_temp[current_mol_type].push_back(line_info);
                    TopParserUtilities::debug_print("Found CMAP entry for ", current_mol_type, ": ", 
                             trimmed, "\n");
                }
            }
        }
        else if (current_section == "molecules") {
            result.molecules_lines.push_back(line_info);
            auto tokens = TopParserUtilities::split(trimmed);
            if (tokens.size() >= 2) {
                TopParserUtilities::debug_print("Found molecule in [ molecules ]: ", tokens[0], 
                         " count=", tokens[1], "\n");
            }
        }
    }

    // Process molecules section
    process_molecules_section(result.molecules_lines, result.molecule_types_order,
                            result.molecule_order, result.used_molecule_types);

    return result;
}

void TopParserSectionProcessor::process_moleculetype_section(const std::string& trimmed, 
                                                           const LineInfo& line_info,
                                                           std::string& current_mol_type,
                                                           std::vector<std::string>& molecule_types_order,
                                                           std::map<std::string, std::vector<LineInfo>>& molecule_atoms_temp,
                                                           std::map<std::string, std::vector<LineInfo>>& molecule_bonds_temp,
                                                           std::map<std::string, std::vector<LineInfo>>& molecule_angles_temp,
                                                           std::map<std::string, std::vector<LineInfo>>& molecule_dihedrals_temp,
                                                           std::map<std::string, std::vector<LineInfo>>& molecule_impropers_temp) {
    (void)line_info; // Parameter intentionally unused currently
    auto tokens = TopParserUtilities::split(trimmed);
    if (tokens.size() >= 2) {
        current_mol_type = tokens[0];
        int current_molecule_nrexcl = std::stoi(tokens[1]);
        
        // If this molecule type already exists, we'll override it
        auto it = std::find(molecule_types_order.begin(), molecule_types_order.end(), 
                          current_mol_type);
        if (it != molecule_types_order.end()) {
            molecule_types_order.erase(it);
            TopParserUtilities::debug_print("Warning: Overriding previous definition of molecule type ", 
                     current_mol_type, "\n");
        }
        molecule_types_order.push_back(current_mol_type);
        
        // Clear all previous definitions for this molecule type
        molecule_atoms_temp[current_mol_type].clear();
        molecule_bonds_temp[current_mol_type].clear();
        molecule_angles_temp[current_mol_type].clear();
        molecule_dihedrals_temp[current_mol_type].clear();
        molecule_impropers_temp[current_mol_type].clear();
        
        TopParserUtilities::debug_print("Processing moleculetype: ", current_mol_type, 
                 "(nrexcl=", current_molecule_nrexcl, ")\n");
    }
}

void TopParserSectionProcessor::process_molecules_section(const std::vector<LineInfo>& molecules_lines,
                                                        const std::vector<std::string>& molecule_types_order,
                                                        std::vector<std::pair<std::string, int>>& molecule_order,
                                                        std::set<std::string>& used_molecule_types) {
    TopParserUtilities::debug_print("\n=== Processing [ molecules ] section ===\n");
    
    if (!molecules_lines.empty()) {
        for (const auto& line_info : molecules_lines) {
            std::string line_no_comment = TopParserUtilities::remove_comment(line_info.content);
            std::string trimmed = TopParserUtilities::trim(line_no_comment);
            if (trimmed.empty()) continue;

            auto tokens = TopParserUtilities::split(trimmed);
            if (tokens.size() >= 2) {
                std::string mol_type = tokens[0];
                int count = std::stoi(tokens[1]);
                molecule_order.push_back(std::make_pair(mol_type, count));
                used_molecule_types.insert(mol_type);
                TopParserUtilities::debug_print("Will add ", count, " copies of molecule ", mol_type, "\n");
            }
        }
    } else if (molecule_types_order.size() == 1) {
        // Special case: if we have exactly one molecule type and no [ molecules ] section,
        // assume we want to add one instance of this molecule (useful for single-molecule itp files)
        std::string mol_type = molecule_types_order[0];
        molecule_order.push_back(std::make_pair(mol_type, 1));
        used_molecule_types.insert(mol_type);
        TopParserUtilities::debug_print("Single molecule type file detected: will add 1 copy of ", mol_type, "\n");
    }
}

} // namespace io
} // namespace pygcmc