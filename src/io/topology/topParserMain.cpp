// src/io/topology/topParserMain.cpp

#include "topParserMain.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <filesystem>
#include <stdexcept>

namespace pygcmc {
namespace io {

model::Topology TOPParser::parse_file(const std::string& filename) {
    model::Topology topology;
    TOPParser parser;
    if (!parser.parse_to_topology(filename, topology)) {
        throw std::runtime_error("Failed to parse topology file: " + filename);
    }
    return topology;
}

model::Topology TOPParser::parse_string(const std::string& top_str) {
    // Handle empty string case - first trim whitespace
    std::string trimmed_str = top_str;
    TopParserUtilities::trim(trimmed_str);
    if (trimmed_str.empty()) {
        return model::Topology();
    }
    
    // Create a temporary file to write the string to
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path();
    std::filesystem::path temp_file = temp_dir / "temp_topology.top";
    
    // Write the string to the temporary file
    std::ofstream out(temp_file);
    if (!out) {
        throw std::runtime_error("Failed to create temporary file for topology parsing");
    }
    out << top_str;  // Write original string to preserve formatting
    out.close();
    
    try {
        // Parse the temporary file
        model::Topology topology = parse_file(temp_file.string());
        
        // Clean up
        std::filesystem::remove(temp_file);
        
        return topology;
    } catch (const std::exception& e) {
        // Clean up on error
        std::filesystem::remove(temp_file);
        throw;
    }
}

bool TOPParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    // Track molecule definitions in order of appearance
    std::vector<std::string> molecule_types_order;
    std::map<std::string, std::vector<LineInfo>> molecule_atoms_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_bonds_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_angles_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_dihedrals_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_impropers_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_cmaps_temp;
    std::vector<LineInfo> molecules_lines;

    TopParserUtilities::debug_print("\n=== Starting topology parsing of ", filename, " ===\n");

    // Clear any previous state
    molecule_atoms_temp.clear();
    molecule_bonds_temp.clear();
    molecule_angles_temp.clear();
    molecule_dihedrals_temp.clear();
    molecule_impropers_temp.clear();
    molecule_cmaps_temp.clear();
    molecule_types_order.clear();
    molecule_order_.clear();
    processed_files_.clear();
    molecules_lines.clear();
    current_molecule_type_.clear();
    current_molecule_nrexcl_ = 0;

    // Track if we found any valid sections
    bool found_valid_section = false;

    // Collect all lines with preprocessor handling
    std::vector<LineInfo> all_lines;
    PreprocessorState pp_state;
    // Check if force field defines are already present in files
    if (!TopParserPreprocessor::collect_all_lines(filename, all_lines, pp_state, true, processed_files_)) {
        return false;
    }
    
    if (all_lines.empty()) {
        TopParserUtilities::debug_print("Error: No valid content found in topology file\n");
        return false;
    }

    // 4) Parse sections using state machine
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
                found_valid_section = true;  // Found a valid section
                TopParserUtilities::debug_print("\n=== Found [ moleculetype ] section at ", line_info.source_file, 
                         ":" , line_info.line_number, " ===\n");
            }
            else if (section_name == "molecules") {
                inside_molecule = false;
                current_mol_type.clear();
                current_section = "molecules";
                found_valid_section = true;  // Found a valid section
                TopParserUtilities::debug_print("\n=== Found [ molecules ] section at ", line_info.source_file, 
                         ":" , line_info.line_number, " ===\n");
            }
            else if (inside_molecule) {
                current_section = section_name;
                found_valid_section = true;  // Found a valid section
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
                auto tokens = TopParserUtilities::split(trimmed);
                if (tokens.size() >= 2) {
                    current_mol_type = tokens[0];
                    current_molecule_nrexcl_ = std::stoi(tokens[1]);
                    
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
                             "(nrexcl=", current_molecule_nrexcl_, ")\n");
                }
            }
            else if (!current_mol_type.empty()) {
                if (current_section == "atoms") {
                    molecule_atoms_temp[current_mol_type].push_back(line_info);
                    auto tokens = TopParserUtilities::split(trimmed);
                    if (tokens.size() >= 4) {  // Assuming at least: nr type resnr resname
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
                    molecule_bonds_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "angles") {
                    molecule_angles_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "dihedrals") {
                    molecule_dihedrals_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "impropers") {
                    molecule_impropers_temp[current_mol_type].push_back(line_info);
                }
                else if (current_section == "cmap") {
                    molecule_cmaps_temp[current_mol_type].push_back(line_info);
                    TopParserUtilities::debug_print("Found CMAP entry for ", current_mol_type, ": ", 
                             trimmed, "\n");
                }
            }
        }
        else if (current_section == "molecules") {
            molecules_lines.push_back(line_info);
            auto tokens = TopParserUtilities::split(trimmed);
            if (tokens.size() >= 2) {
                TopParserUtilities::debug_print("Found molecule in [ molecules ]: ", tokens[0], 
                         " count=", tokens[1], "\n");
            }
        }
    }

    // Return false if no valid sections were found
    if (!found_valid_section) {
        TopParserUtilities::debug_print("Error: No valid topology sections found in file\n");
        return false;
    }

    // After collecting all definitions, process the [ molecules ] section first
    TopParserUtilities::debug_print("\n=== Processing [ molecules ] section ===\n");
    std::set<std::string> used_molecule_types;
    if (!molecules_lines.empty()) {
        for (const auto& line_info : molecules_lines) {
            std::string line_no_comment = TopParserUtilities::remove_comment(line_info.content);
            std::string trimmed = TopParserUtilities::trim(line_no_comment);
            if (trimmed.empty()) continue;

            auto tokens = TopParserUtilities::split(trimmed);
            if (tokens.size() >= 2) {
                std::string mol_type = tokens[0];
                int count = std::stoi(tokens[1]);
                molecule_order_.push_back(std::make_pair(mol_type, count));
                used_molecule_types.insert(mol_type);
                TopParserUtilities::debug_print("Will add ", count, " copies of molecule ", mol_type, "\n");
            }
        }
    } else if (molecule_types_order.size() == 1) {
        // Special case: if we have exactly one molecule type and no [ molecules ] section,
        // assume we want to add one instance of this molecule (useful for single-molecule itp files)
        std::string mol_type = molecule_types_order[0];
        molecule_order_.push_back(std::make_pair(mol_type, 1));
        used_molecule_types.insert(mol_type);
        TopParserUtilities::debug_print("Single molecule type file detected: will add 1 copy of ", mol_type, "\n");
    }

    // Now we know which molecule types are actually used, we can write them to topology
    TopParserUtilities::debug_print("\n=== Counting atoms per molecule type ===\n");
    int atom_offset = 0;
    int mol_index = 0;
    std::map<std::string, int> atoms_per_molecule;

    // First count atoms per molecule type, but only for molecules we'll actually use
    for (const auto& mol_type : used_molecule_types) {
        if (molecule_atoms_temp.find(mol_type) != molecule_atoms_temp.end()) {
            int atom_count = 0;
            TopParserUtilities::debug_print("Checking atoms for molecule type ", mol_type, ":\n");
            for (const auto& line_info : molecule_atoms_temp[mol_type]) {
                auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line_info.content));
                if (tokens.size() >= 8) {
                    atom_count++;
                    TopParserUtilities::debug_print("  Atom ", tokens[0], " (", tokens[1], 
                             ") in residue ", tokens[3], "\n");
                }
            }
            atoms_per_molecule[mol_type] = atom_count;
            TopParserUtilities::debug_print("Molecule ", mol_type, " has ", atom_count, " atoms\n");
        } else {
            TopParserUtilities::debug_print("Warning: No atom definitions found for molecule type ", mol_type, "\n");
        }
    }

    TopParserUtilities::debug_print("\n=== Adding molecules to topology ===\n");
    // Then add molecules in order according to [ molecules ] section
    for (const auto& mol : molecule_order_) {
        const std::string& mol_type = mol.first;
        int count = mol.second;

        // Skip if we don't have the molecule definition
        if (molecule_atoms_temp.find(mol_type) == molecule_atoms_temp.end()) {
            TopParserUtilities::debug_print("Warning: No definition found for molecule type ", mol_type, "\n");
            continue;
        }

        TopParserUtilities::debug_print("Adding ", count, " copies of molecule ", mol_type, "\n");

        // Add 'count' copies of this molecule type
        for (int i = 0; i < count; i++) {
            // Create unique segment name for this molecule instance
            std::string segment_name = mol_type + "_" + std::to_string(mol_index++);
            topology.add_segment(segment_name);
            current_molecule_type_ = segment_name;
            TopParserUtilities::debug_print("  Creating segment ", segment_name, "\n");

            // Add atoms for this molecule instance
            if (!TopParserSections::parse_atoms_section(molecule_atoms_temp[mol_type], topology, segment_name)) {
                const auto& line = molecule_atoms_temp[mol_type].front();
                TopParserUtilities::debug_print("Error parsing atoms for molecule ", mol_type, 
                         " at ", line.source_file, ":", line.line_number, "\n");
                return false;
            }
            TopParserUtilities::debug_print("  Added atoms for segment ", segment_name, "\n");

            // Add bonds for this molecule instance
            if (molecule_bonds_temp.find(mol_type) != molecule_bonds_temp.end()) {
                if (!TopParserSections::parse_bonds_section(molecule_bonds_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_bonds_temp[mol_type].front();
                    TopParserUtilities::debug_print("Error parsing bonds for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                TopParserUtilities::debug_print("  Added bonds for segment ", segment_name, "\n");
            }

            // Add angles for this molecule instance
            if (molecule_angles_temp.find(mol_type) != molecule_angles_temp.end()) {
                if (!TopParserSections::parse_angles_section(molecule_angles_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_angles_temp[mol_type].front();
                    TopParserUtilities::debug_print("Error parsing angles for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                TopParserUtilities::debug_print("  Added angles for segment ", segment_name, "\n");
            }

            // Add dihedrals for this molecule instance
            if (molecule_dihedrals_temp.find(mol_type) != molecule_dihedrals_temp.end()) {
                if (!TopParserSections::parse_dihedrals_section(molecule_dihedrals_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_dihedrals_temp[mol_type].front();
                    TopParserUtilities::debug_print("Error parsing dihedrals for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                TopParserUtilities::debug_print("  Added dihedrals for segment ", segment_name, "\n");
            }

            // Add impropers for this molecule instance
            if (molecule_impropers_temp.find(mol_type) != molecule_impropers_temp.end()) {
                if (!TopParserSections::parse_impropers_section(molecule_impropers_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_impropers_temp[mol_type].front();
                    TopParserUtilities::debug_print("Error parsing impropers for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                TopParserUtilities::debug_print("  Added impropers for segment ", segment_name, "\n");
            }

            // Add CMAPs for this molecule instance
            if (molecule_cmaps_temp.find(mol_type) != molecule_cmaps_temp.end()) {
                if (!TopParserSections::parse_cmaps_section(molecule_cmaps_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_cmaps_temp[mol_type].front();
                    TopParserUtilities::debug_print("Error parsing CMAPs for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                TopParserUtilities::debug_print("  Added CMAPs for segment ", segment_name, "\n");
            }

            // Update atom offset for next molecule instance
            atom_offset += atoms_per_molecule[mol_type];
        }
    }

    TopParserUtilities::debug_print("\n=== Final topology statistics ===\n");
    TopParserUtilities::debug_print("Total atoms: ", topology.get_num_atoms(), "\n");
    TopParserUtilities::debug_print("Total bonds: ", topology.get_num_bonds(), "\n");
    TopParserUtilities::debug_print("Total angles: ", topology.get_num_angles(), "\n");
    TopParserUtilities::debug_print("Total dihedrals: ", topology.get_num_dihedrals(), "\n");
    TopParserUtilities::debug_print("Total impropers: ", topology.get_num_impropers(), "\n");
    TopParserUtilities::debug_print("Total residues: ", topology.get_num_residues(), "\n");
    TopParserUtilities::debug_print("Total segments: ", topology.get_num_segments(), "\n");

    return true;
}

} // namespace io
} // namespace pygcmc