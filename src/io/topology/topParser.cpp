// src/io/topParser.cpp

#include "topParser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <filesystem>
#include <stdexcept>

namespace pygcmc {
namespace io {

// Initialize static member
bool TOPParser::debug_enabled_ = false;

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
    trim(trimmed_str);
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

    debug_print("\n=== Starting topology parsing of ", filename, " ===\n");

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
    if (!collect_all_lines(filename, all_lines, pp_state, true)) {
        return false;
    }
    
    if (all_lines.empty()) {
        debug_print("Error: No valid content found in topology file\n");
        return false;
    }

    // 4) Parse sections using state machine
    std::string current_mol_type;
    bool inside_molecule = false;
    std::string current_section;

    for (const auto& line_info : all_lines) {
        std::string line_no_comment = remove_comment(line_info.content);
        std::string trimmed = trim(line_no_comment);
        if (trimmed.empty()) {
            continue;
        }

        // Check if line starts with '['
        if (trimmed.size() > 1 && trimmed[0] == '[' && trimmed[trimmed.size() - 1] == ']') {
            // Remove the leading '[' and trailing ']' and trim whitespace
            std::string section_name = trimmed.substr(1, trimmed.size() - 2);
            section_name = trim(section_name);

            // Convert to lowercase
            std::transform(section_name.begin(), section_name.end(), 
                         section_name.begin(), ::tolower);

            if (section_name == "moleculetype") {
                inside_molecule = true;
                current_mol_type.clear();
                current_section = "moleculetype";
                found_valid_section = true;  // Found a valid section
                debug_print("\n=== Found [ moleculetype ] section at ", line_info.source_file, 
                         ":" , line_info.line_number, " ===\n");
            }
            else if (section_name == "molecules") {
                inside_molecule = false;
                current_mol_type.clear();
                current_section = "molecules";
                found_valid_section = true;  // Found a valid section
                debug_print("\n=== Found [ molecules ] section at ", line_info.source_file, 
                         ":" , line_info.line_number, " ===\n");
            }
            else if (inside_molecule) {
                current_section = section_name;
                found_valid_section = true;  // Found a valid section
                debug_print("Found [ ", section_name, " ] section for molecule ", 
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
                auto tokens = split(trimmed);
                if (tokens.size() >= 2) {
                    current_mol_type = tokens[0];
                    current_molecule_nrexcl_ = std::stoi(tokens[1]);
                    
                    // If this molecule type already exists, we'll override it
                    auto it = std::find(molecule_types_order.begin(), molecule_types_order.end(), 
                                      current_mol_type);
                    if (it != molecule_types_order.end()) {
                        molecule_types_order.erase(it);
                        debug_print("Warning: Overriding previous definition of molecule type ", 
                                 current_mol_type, "\n");
                    }
                    molecule_types_order.push_back(current_mol_type);
                    
                    // Clear all previous definitions for this molecule type
                    molecule_atoms_temp[current_mol_type].clear();
                    molecule_bonds_temp[current_mol_type].clear();
                    molecule_angles_temp[current_mol_type].clear();
                    molecule_dihedrals_temp[current_mol_type].clear();
                    molecule_impropers_temp[current_mol_type].clear();
                    
                    debug_print("Processing moleculetype: ", current_mol_type, 
                             "(nrexcl=", current_molecule_nrexcl_, ")\n");
                }
            }
            else if (!current_mol_type.empty()) {
                if (current_section == "atoms") {
                    molecule_atoms_temp[current_mol_type].push_back(line_info);
                    auto tokens = split(trimmed);
                    if (tokens.size() >= 4) {  // Assuming at least: nr type resnr resname
                        debug_print("Added atom to ", current_mol_type, ": ", 
                                 tokens[0], " ", tokens[1], " (residue ", 
                                 tokens[3], ")\n");
                    }
                }
                else if (current_section == "settles") {
                    debug_print("Found settles for ", current_mol_type, ": ", 
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
                    debug_print("Found CMAP entry for ", current_mol_type, ": ", 
                             trimmed, "\n");
                }
            }
        }
        else if (current_section == "molecules") {
            molecules_lines.push_back(line_info);
            auto tokens = split(trimmed);
            if (tokens.size() >= 2) {
                debug_print("Found molecule in [ molecules ]: ", tokens[0], 
                         " count=", tokens[1], "\n");
            }
        }
    }

    // Return false if no valid sections were found
    if (!found_valid_section) {
        debug_print("Error: No valid topology sections found in file\n");
        return false;
    }

    // After collecting all definitions, process the [ molecules ] section first
    debug_print("\n=== Processing [ molecules ] section ===\n");
    std::set<std::string> used_molecule_types;
    if (!molecules_lines.empty()) {
        for (const auto& line_info : molecules_lines) {
            std::string line_no_comment = remove_comment(line_info.content);
            std::string trimmed = trim(line_no_comment);
            if (trimmed.empty()) continue;

            auto tokens = split(trimmed);
            if (tokens.size() >= 2) {
                std::string mol_type = tokens[0];
                int count = std::stoi(tokens[1]);
                molecule_order_.push_back(std::make_pair(mol_type, count));
                used_molecule_types.insert(mol_type);
                debug_print("Will add ", count, " copies of molecule ", mol_type, "\n");
            }
        }
    } else if (molecule_types_order.size() == 1) {
        // Special case: if we have exactly one molecule type and no [ molecules ] section,
        // assume we want to add one instance of this molecule (useful for single-molecule itp files)
        std::string mol_type = molecule_types_order[0];
        molecule_order_.push_back(std::make_pair(mol_type, 1));
        used_molecule_types.insert(mol_type);
        debug_print("Single molecule type file detected: will add 1 copy of ", mol_type, "\n");
    }

    // Now we know which molecule types are actually used, we can write them to topology
    debug_print("\n=== Counting atoms per molecule type ===\n");
    int atom_offset = 0;
    int mol_index = 0;
    std::map<std::string, int> atoms_per_molecule;

    // First count atoms per molecule type, but only for molecules we'll actually use
    for (const auto& mol_type : used_molecule_types) {
        if (molecule_atoms_temp.find(mol_type) != molecule_atoms_temp.end()) {
            int atom_count = 0;
            debug_print("Checking atoms for molecule type ", mol_type, ":\n");
            for (const auto& line_info : molecule_atoms_temp[mol_type]) {
                auto tokens = split(remove_comment(line_info.content));
                if (tokens.size() >= 8) {
                    atom_count++;
                    debug_print("  Atom ", tokens[0], " (", tokens[1], 
                             ") in residue ", tokens[3], "\n");
                }
            }
            atoms_per_molecule[mol_type] = atom_count;
            debug_print("Molecule ", mol_type, " has ", atom_count, " atoms\n");
        } else {
            debug_print("Warning: No atom definitions found for molecule type ", mol_type, "\n");
        }
    }

    debug_print("\n=== Adding molecules to topology ===\n");
    // Then add molecules in order according to [ molecules ] section
    for (const auto& mol : molecule_order_) {
        const std::string& mol_type = mol.first;
        int count = mol.second;

        // Skip if we don't have the molecule definition
        if (molecule_atoms_temp.find(mol_type) == molecule_atoms_temp.end()) {
            debug_print("Warning: No definition found for molecule type ", mol_type, "\n");
            continue;
        }

        debug_print("Adding ", count, " copies of molecule ", mol_type, "\n");

        // Add 'count' copies of this molecule type
        for (int i = 0; i < count; i++) {
            // Create unique segment name for this molecule instance
            std::string segment_name = mol_type + "_" + std::to_string(mol_index++);
            topology.add_segment(segment_name);
            current_molecule_type_ = segment_name;
            debug_print("  Creating segment ", segment_name, "\n");

            // Add atoms for this molecule instance
            if (!parse_atoms_section(molecule_atoms_temp[mol_type], topology)) {
                const auto& line = molecule_atoms_temp[mol_type].front();
                debug_print("Error parsing atoms for molecule ", mol_type, 
                         " at ", line.source_file, ":", line.line_number, "\n");
                return false;
            }
            debug_print("  Added atoms for segment ", segment_name, "\n");

            // Add bonds for this molecule instance
            if (molecule_bonds_temp.find(mol_type) != molecule_bonds_temp.end()) {
                if (!parse_bonds_section(molecule_bonds_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_bonds_temp[mol_type].front();
                    debug_print("Error parsing bonds for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                debug_print("  Added bonds for segment ", segment_name, "\n");
            }

            // Add angles for this molecule instance
            if (molecule_angles_temp.find(mol_type) != molecule_angles_temp.end()) {
                if (!parse_angles_section(molecule_angles_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_angles_temp[mol_type].front();
                    debug_print("Error parsing angles for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                debug_print("  Added angles for segment ", segment_name, "\n");
            }

            // Add dihedrals for this molecule instance
            if (molecule_dihedrals_temp.find(mol_type) != molecule_dihedrals_temp.end()) {
                if (!parse_dihedrals_section(molecule_dihedrals_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_dihedrals_temp[mol_type].front();
                    debug_print("Error parsing dihedrals for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                debug_print("  Added dihedrals for segment ", segment_name, "\n");
            }

            // Add impropers for this molecule instance
            if (molecule_impropers_temp.find(mol_type) != molecule_impropers_temp.end()) {
                if (!parse_impropers_section(molecule_impropers_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_impropers_temp[mol_type].front();
                    debug_print("Error parsing impropers for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                debug_print("  Added impropers for segment ", segment_name, "\n");
            }

            // Add CMAPs for this molecule instance
            if (molecule_cmaps_temp.find(mol_type) != molecule_cmaps_temp.end()) {
                if (!parse_cmaps_section(molecule_cmaps_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_cmaps_temp[mol_type].front();
                    debug_print("Error parsing CMAPs for molecule ", mol_type, 
                             " at ", line.source_file, ":", line.line_number, "\n");
                    return false;
                }
                debug_print("  Added CMAPs for segment ", segment_name, "\n");
            }

            // Update atom offset for next molecule instance
            atom_offset += atoms_per_molecule[mol_type];
        }
    }

    debug_print("\n=== Final topology statistics ===\n");
    debug_print("Total atoms: ", topology.get_num_atoms(), "\n");
    debug_print("Total bonds: ", topology.get_num_bonds(), "\n");
    debug_print("Total angles: ", topology.get_num_angles(), "\n");
    debug_print("Total dihedrals: ", topology.get_num_dihedrals(), "\n");
    debug_print("Total impropers: ", topology.get_num_impropers(), "\n");
    debug_print("Total residues: ", topology.get_num_residues(), "\n");
    debug_print("Total segments: ", topology.get_num_segments(), "\n");

    return true;
}

bool TOPParser::collect_all_lines(const std::string& filename, std::vector<LineInfo>& all_lines,
                                PreprocessorState& pp_state, bool is_main_file) {
    // Check if we've already processed this file
    if (processed_files_.find(filename) != processed_files_.end()) {
        return true;
    }
    processed_files_.insert(filename);

    std::ifstream file(filename);
    if (!file.is_open()) {
        if (is_main_file) {
            debug_print("Error: Cannot open main topology file '", filename, "'\n");
            return false;
        }
        debug_print("Warning: Cannot open included file '", filename, "'\n");
        return true;  // Continue for included files
    }

    std::string line;
    int line_number = 0;

    while (std::getline(file, line)) {
        line_number++;
        trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == ';') {
            continue;
        }

        // Handle preprocessor directives
        if (line[0] == '#') {
            if (!process_preprocessor_line(line, filename, all_lines, pp_state, line_number)) {
                return false;
            }
            continue;
        }

        // Skip lines if inside a false #ifdef/#ifndef block
        if (pp_state.skip_section) {
            debug_print("Skipping line due to preprocessor: ", line, "\n");
            continue;
        }

        // Add normal content line with source information
        all_lines.emplace_back(line, filename, line_number);
    }

    return true;
}

bool TOPParser::process_preprocessor_line(const std::string& line, const std::string& parent_file,
                                        std::vector<LineInfo>& all_lines, PreprocessorState& pp_state,
                                        int line_number) {
    std::istringstream iss(line);
    std::string directive;
    iss >> directive;

    if (directive == "#include") {
        // Only process include if we're not in a false ifdef block
        if (!pp_state.should_skip()) {
            std::string include_path;
            std::getline(iss, include_path);
            trim(include_path);
            
            // Remove quotes if present
            if (include_path.front() == '"' && include_path.back() == '"') {
                include_path = include_path.substr(1, include_path.length() - 2);
            }

            std::string resolved_path = resolve_include_path(include_path, parent_file);
            if (!resolved_path.empty()) {
                if (!collect_all_lines(resolved_path, all_lines, pp_state, false)) {
                    return false;
                }
            } else {
                debug_print("Warning: Include file not found: ", include_path, 
                         " (referenced from ", parent_file, ":", line_number, ")\n");
            }
        }
    }
    else if (directive == "#ifdef" || directive == "#ifndef") {
        std::string macro_name;
        iss >> macro_name;
        bool is_defined = (pp_state.defines.find(macro_name) != pp_state.defines.end());
        
        // If we're already in a skipped section, push false to maintain nesting
        if (pp_state.should_skip()) {
            pp_state.ifdef_stack.push_back(false);
            pp_state.else_encountered.push_back(false);
        } else {
            bool condition_met = (directive == "#ifdef") ? is_defined : !is_defined;
            pp_state.ifdef_stack.push_back(condition_met);
            pp_state.else_encountered.push_back(false);
        }
        
        pp_state.skip_section = pp_state.should_skip();
        
        debug_print("Processing ", directive, " ", macro_name, 
                 ": defined=", is_defined, ", skip=", pp_state.skip_section, 
                 ", stack_size=", pp_state.ifdef_stack.size(), "\n");
    }
    else if (directive == "#else") {
        if (pp_state.ifdef_stack.empty()) {
            debug_print("Warning: Unmatched #else at ", parent_file, ":", line_number, "\n");
            return false;
        }
        
        // Only flip the condition if we haven't seen an #else at this level yet
        // AND we're not in an outer skipped section
        if (!pp_state.else_encountered.back()) {
            bool outer_skip = false;
            // Check if any outer level is false (skipped)
            for (size_t i = 0; i < pp_state.ifdef_stack.size() - 1; ++i) {
                if (!pp_state.ifdef_stack[i]) {
                    outer_skip = true;
                    break;
                }
            }
            
            if (!outer_skip) {
                pp_state.ifdef_stack.back() = !pp_state.ifdef_stack.back();
            }
            pp_state.else_encountered.back() = true;
            pp_state.skip_section = pp_state.should_skip();
            
            debug_print("Processing #else: skip=", pp_state.skip_section, 
                     ", stack_size=", pp_state.ifdef_stack.size(), 
                     ", outer_skip=", outer_skip, "\n");
        } else {
            debug_print("Warning: Multiple #else directives at the same nesting level at ",
                     parent_file, ":", line_number, "\n");
        }
    }
    else if (directive == "#endif") {
        if (pp_state.ifdef_stack.empty()) {
            debug_print("Warning: Unmatched #endif at ", parent_file, ":", line_number, "\n");
            return false;
        }
        
        pp_state.ifdef_stack.pop_back();
        pp_state.else_encountered.pop_back();
        pp_state.skip_section = pp_state.should_skip();
        
        debug_print("Processing #endif: skip=", pp_state.skip_section, 
                 ", stack_size=", pp_state.ifdef_stack.size(), "\n");
    }
    else if (directive == "#define") {
        if (!pp_state.should_skip()) {
            std::string macro_name;
            iss >> macro_name;
            std::string macro_value;
            std::getline(iss, macro_value);
            trim(macro_value);
            pp_state.defines[macro_name] = macro_value;
            debug_print("Defined macro: ", macro_name, " = ", macro_value, "\n");
        }
    }
    else if (directive == "#undef") {
        if (!pp_state.should_skip()) {
            std::string macro_name;
            iss >> macro_name;
            pp_state.defines.erase(macro_name);
            debug_print("Undefined macro: ", macro_name, "\n");
        }
    }
    
    return true;
}

void TOPParser::parse_sections(const std::vector<LineInfo>& all_lines,
                             std::map<std::string, std::vector<LineInfo>>& sections) {
    std::string current_section;
    for (const auto& line_info : all_lines) {
        const std::string& line = line_info.content;
        
        // Check for section header
        if (line[0] == '[') {
            std::string tmp = line;
            trim(tmp);
            // Remove brackets and trim again
            if (tmp.front() == '[') tmp.erase(tmp.begin());
            if (!tmp.empty() && tmp.back() == ']') tmp.pop_back();
            trim(tmp);
            current_section = tmp;
            
            // Ensure section exists in map
            if (sections.find(current_section) == sections.end()) {
                sections[current_section] = std::vector<LineInfo>();
            }
            continue;
        }

        // Add line to current section if we're in one
        if (!current_section.empty()) {
            sections[current_section].push_back(line_info);
        }
    }
}

std::string TOPParser::remove_comment(const std::string& line) {
    size_t comment_pos = line.find(';');
    if (comment_pos != std::string::npos) {
        return line.substr(0, comment_pos);
    }
    return line;
}

bool TOPParser::parse_moleculetype_section(const std::vector<LineInfo>& lines, [[maybe_unused]] model::Topology& topology) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
        if (tokens.size() < 2) continue;

        try {
            current_molecule_type_ = tokens[0];
            current_molecule_nrexcl_ = std::stoi(tokens[1]);
            return true;
        } catch (const std::exception& e) {
            debug_print("Error parsing moleculetype line: ", line, "\n");
            return false;
        }
    }
    return true;
}

bool TOPParser::parse_atoms_section(const std::vector<LineInfo>& lines, model::Topology& topology) {
    // Get current residue number offset from segment name
    int segment_index = 0;
    size_t underscore_pos = current_molecule_type_.find_last_of('_');
    if (underscore_pos != std::string::npos) {
        try {
            segment_index = std::stoi(current_molecule_type_.substr(underscore_pos + 1));
        } catch (...) {
            segment_index = 0;
        }
    }

    // Find the maximum residue number in this molecule definition
    int max_resnum = 0;
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
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
        auto tokens = split(remove_comment(line));
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
                                              : default_mass_for_atom_type(atom_type);

            // Adjust residue number based on segment index
            int adjusted_resnum = residue_number + (segment_index * max_resnum);

            // Add atom to topology using current_molecule_type_ as segment
            topology.add_atom(
                atom_name,
                atom_type,
                charge,
                mass,
                residue_name,
                adjusted_resnum,
                current_molecule_type_
            );
        } catch (const std::exception& e) {
            debug_print("Error parsing atom line: ", line, "\n");
            return false;
        }
    }
    return true;
}

bool TOPParser::parse_bonds_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
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
            debug_print("Warning: Error parsing bond line: ", line, " at ", 
                     line_info.source_file, ":", line_info.line_number, 
                     " - ", e.what(), "\n");
            continue;  // Continue with next line instead of failing
        }
    }
    return true;  // Return true if we processed all lines (even with some warnings)
}

bool TOPParser::parse_angles_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
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
            debug_print("Warning: Error parsing angle line: ", line, " at ",
                     line_info.source_file, ":", line_info.line_number, 
                     " - ", e.what(), "\n");
            continue;  // Continue with next line instead of failing
        }
    }
    return true;
}

bool TOPParser::parse_dihedrals_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    debug_print("Parsing ", lines.size(), " dihedral lines\n");
    int proper_count = 0;
    int improper_count = 0;

    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
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
            debug_print("Error parsing dihedral line: ", line, " - ", e.what(), "\n");
            return false;
        }
    }
    
    debug_print("Added ", proper_count, " proper dihedrals and ", 
              improper_count, " improper dihedrals\n");
    return true;
}

bool TOPParser::parse_impropers_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
        if (tokens.size() < 4) continue;  // Need at least ai, aj, ak, al

        try {
            int atom1 = std::stoi(tokens[0]) - 1 + atom_offset;
            int atom2 = std::stoi(tokens[1]) - 1 + atom_offset;
            int atom3 = std::stoi(tokens[2]) - 1 + atom_offset;
            int atom4 = std::stoi(tokens[3]) - 1 + atom_offset;
            
            // Add improper to topology
            topology.add_improper(atom1, atom2, atom3, atom4);
        } catch (const std::exception& e) {
            debug_print("Error parsing improper line: ", line, "\n");
            return false;
        }
    }
    return true;
}

bool TOPParser::parse_molecules_section(const std::vector<LineInfo>& lines, [[maybe_unused]] model::Topology& topology) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
        if (tokens.size() < 2) continue;

        try {
            std::string mol_type = tokens[0];
            int count = std::stoi(tokens[1]);
            
            // Store molecule type and count
            molecule_order_.push_back(std::make_pair(mol_type, count));
        } catch (const std::exception& e) {
            debug_print("Error parsing molecules line: ", line, "\n");
            return false;
        }
    }
    return true;
}

std::string TOPParser::trim(std::string& str) {
    // Trim leading spaces
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    
    // Trim trailing spaces
    str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), str.end());
    
    return str;
}

std::vector<std::string> TOPParser::split(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    
    while (iss >> token) {
        if (token[0] == ';') break;  // Stop at comments
        tokens.push_back(token);
    }
    
    return tokens;
}

std::string TOPParser::resolve_include_path(const std::string& include_path, const std::string& parent_file) {
    namespace fs = std::filesystem;
    
    // Convert parent path to absolute and get its directory
    fs::path parent_path = fs::absolute(parent_file);
    fs::path parent_dir = parent_path.parent_path();
    
    // Debug output
    debug_print("Resolving include path: ", include_path, "\n",
             "Parent file: ", parent_file, "\n",
             "Parent dir: ", parent_dir.string(), "\n");
    
    // Simply combine parent directory with include path
    fs::path resolved = parent_dir / include_path;
    if (fs::exists(resolved)) {
        debug_print("Found include file at: ", resolved.string(), "\n");
        return resolved.string();
    }
    
    // If not found, provide error message
    debug_print("Warning: Include file not found: ", include_path, "\n",
             "Tried path: ", resolved.string(), "\n");
    
    return "";
}

bool TOPParser::parse_cmaps_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
        if (tokens.size() < 6) {  // Need 5 atoms + function type
            debug_print("Warning: Skipping CMAP line with insufficient tokens: ", line, "\n");
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
            
            debug_print("Added CMAP between atoms: ", cmap_atoms[0], " ", cmap_atoms[1], " ", cmap_atoms[2], " ", cmap_atoms[3], " ", cmap_atoms[4], " (function type ", function_type, ")\n");
        } catch (const std::exception& e) {
            debug_print("Error parsing CMAP line: ", line, " at ", 
                     line_info.source_file, ":", line_info.line_number, 
                     " - ", e.what(), "\n");
            return false;
        }
    }
    return true;
}

// 新增: 根据原子类型推断默认质量
double TOPParser::default_mass_for_atom_type(const std::string& atom_type) {
    auto starts_with = [](const std::string& s, const std::string& prefix) {
        return s.rfind(prefix, 0) == 0;
    };

    if (atom_type.empty()) {
        return 1.0; // fallback
    }

    switch (atom_type[0]) {
        case 'H':
            return 1.008; // Hydrogen
        case 'C':
            if (starts_with(atom_type, "CT") || starts_with(atom_type, "CA")) {
                return 12.011;
            }
            return 12.011; // Carbon
        case 'N':
            if (starts_with(atom_type, "NH")) {
                return 14.007;
            }
            return 14.007; // Nitrogen
        case 'O':
            if (starts_with(atom_type, "OT") || starts_with(atom_type, "OH")) {
                return 15.999;
            }
            return 15.999; // Oxygen
        case 'S':
            return 32.065; // Sulfur
        case 'P':
            return 30.974; // Phosphorus
        case 'L':
            return 0.0;    // Lone pair virtual sites
        case 'D':
            return 0.4;    // Drude oscillator particles
        default:
            return 1.0;    // Generic fallback for unknown types
    }
}

} // namespace io
} // namespace pygcmc





