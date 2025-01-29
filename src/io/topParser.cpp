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

bool TopParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    // Track molecule definitions in order of appearance
    std::vector<std::string> molecule_types_order;
    std::map<std::string, std::vector<LineInfo>> molecule_atoms_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_bonds_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_angles_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_dihedrals_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_impropers_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_cmaps_temp;
    std::vector<LineInfo> molecules_lines;

    std::cerr << "\n=== Starting topology parsing of " << filename << " ===" << std::endl;

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
    if (!collect_all_lines(filename, all_lines, pp_state, true)) {
        return false;
    }
    
    if (all_lines.empty()) {
        std::cerr << "Error: No valid content found in topology file" << std::endl;
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
                std::cerr << "\n=== Found [ moleculetype ] section at " << line_info.source_file 
                         << ":" << line_info.line_number << " ===" << std::endl;
            }
            else if (section_name == "molecules") {
                inside_molecule = false;
                current_mol_type.clear();
                current_section = "molecules";
                found_valid_section = true;  // Found a valid section
                std::cerr << "\n=== Found [ molecules ] section at " << line_info.source_file 
                         << ":" << line_info.line_number << " ===" << std::endl;
            }
            else if (inside_molecule) {
                current_section = section_name;
                found_valid_section = true;  // Found a valid section
                std::cerr << "Found [ " << section_name << " ] section for molecule " 
                         << current_mol_type << std::endl;
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
                        std::cerr << "Warning: Overriding previous definition of molecule type " 
                                 << current_mol_type << std::endl;
                    }
                    molecule_types_order.push_back(current_mol_type);
                    
                    // Clear all previous definitions for this molecule type
                    molecule_atoms_temp[current_mol_type].clear();
                    molecule_bonds_temp[current_mol_type].clear();
                    molecule_angles_temp[current_mol_type].clear();
                    molecule_dihedrals_temp[current_mol_type].clear();
                    molecule_impropers_temp[current_mol_type].clear();
                    
                    std::cerr << "Processing moleculetype: " << current_mol_type 
                             << " (nrexcl=" << current_molecule_nrexcl_ << ")" << std::endl;
                }
            }
            else if (!current_mol_type.empty()) {
                if (current_section == "atoms") {
                    molecule_atoms_temp[current_mol_type].push_back(line_info);
                    auto tokens = split(trimmed);
                    if (tokens.size() >= 4) {  // Assuming at least: nr type resnr resname
                        std::cerr << "Added atom to " << current_mol_type << ": " 
                                 << tokens[0] << " " << tokens[1] << " (residue " 
                                 << tokens[3] << ")" << std::endl;
                    }
                }
                else if (current_section == "settles") {
                    std::cerr << "Found settles for " << current_mol_type << ": " 
                             << trimmed << std::endl;
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
                    std::cerr << "Found CMAP entry for " << current_mol_type << ": " 
                             << trimmed << std::endl;
                }
            }
        }
        else if (current_section == "molecules") {
            molecules_lines.push_back(line_info);
            auto tokens = split(trimmed);
            if (tokens.size() >= 2) {
                std::cerr << "Found molecule in [ molecules ]: " << tokens[0] 
                         << " count=" << tokens[1] << std::endl;
            }
        }
    }

    // Return false if no valid sections were found
    if (!found_valid_section) {
        std::cerr << "Error: No valid topology sections found in file" << std::endl;
        return false;
    }

    // After collecting all definitions, process the [ molecules ] section first
    std::cerr << "\n=== Processing [ molecules ] section ===" << std::endl;
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
                std::cerr << "Will add " << count << " copies of molecule " << mol_type << std::endl;
            }
        }
    } else if (molecule_types_order.size() == 1) {
        // Special case: if we have exactly one molecule type and no [ molecules ] section,
        // assume we want to add one instance of this molecule (useful for single-molecule itp files)
        std::string mol_type = molecule_types_order[0];
        molecule_order_.push_back(std::make_pair(mol_type, 1));
        used_molecule_types.insert(mol_type);
        std::cerr << "Single molecule type file detected: will add 1 copy of " << mol_type << std::endl;
    }

    // Now we know which molecule types are actually used, we can write them to topology
    std::cerr << "\n=== Counting atoms per molecule type ===" << std::endl;
    int atom_offset = 0;
    int mol_index = 0;
    std::map<std::string, int> atoms_per_molecule;

    // First count atoms per molecule type, but only for molecules we'll actually use
    for (const auto& mol_type : used_molecule_types) {
        if (molecule_atoms_temp.find(mol_type) != molecule_atoms_temp.end()) {
            int atom_count = 0;
            std::cerr << "Checking atoms for molecule type " << mol_type << ":" << std::endl;
            for (const auto& line_info : molecule_atoms_temp[mol_type]) {
                auto tokens = split(remove_comment(line_info.content));
                if (tokens.size() >= 8) {
                    atom_count++;
                    std::cerr << "  Atom " << tokens[0] << " (" << tokens[1] 
                             << ") in residue " << tokens[3] << std::endl;
                }
            }
            atoms_per_molecule[mol_type] = atom_count;
            std::cerr << "Molecule " << mol_type << " has " << atom_count << " atoms" << std::endl;
        } else {
            std::cerr << "Warning: No atom definitions found for molecule type " << mol_type << std::endl;
        }
    }

    std::cerr << "\n=== Adding molecules to topology ===" << std::endl;
    // Then add molecules in order according to [ molecules ] section
    for (const auto& mol : molecule_order_) {
        const std::string& mol_type = mol.first;
        int count = mol.second;

        // Skip if we don't have the molecule definition
        if (molecule_atoms_temp.find(mol_type) == molecule_atoms_temp.end()) {
            std::cerr << "Warning: No definition found for molecule type " << mol_type << std::endl;
            continue;
        }

        std::cerr << "Adding " << count << " copies of molecule " << mol_type << std::endl;

        // Add 'count' copies of this molecule type
        for (int i = 0; i < count; i++) {
            // Create unique segment name for this molecule instance
            std::string segment_name = mol_type + "_" + std::to_string(mol_index++);
            topology.add_segment(segment_name);
            current_molecule_type_ = segment_name;
            std::cerr << "  Creating segment " << segment_name << std::endl;

            // Add atoms for this molecule instance
            if (!parse_atoms_section(molecule_atoms_temp[mol_type], topology)) {
                const auto& line = molecule_atoms_temp[mol_type].front();
                std::cerr << "Error parsing atoms for molecule " << mol_type 
                         << " at " << line.source_file << ":" << line.line_number << std::endl;
                return false;
            }
            std::cerr << "  Added atoms for segment " << segment_name << std::endl;

            // Add bonds for this molecule instance
            if (molecule_bonds_temp.find(mol_type) != molecule_bonds_temp.end()) {
                if (!parse_bonds_section(molecule_bonds_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_bonds_temp[mol_type].front();
                    std::cerr << "Error parsing bonds for molecule " << mol_type 
                             << " at " << line.source_file << ":" << line.line_number << std::endl;
                    return false;
                }
                std::cerr << "  Added bonds for segment " << segment_name << std::endl;
            }

            // Add angles for this molecule instance
            if (molecule_angles_temp.find(mol_type) != molecule_angles_temp.end()) {
                if (!parse_angles_section(molecule_angles_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_angles_temp[mol_type].front();
                    std::cerr << "Error parsing angles for molecule " << mol_type 
                             << " at " << line.source_file << ":" << line.line_number << std::endl;
                    return false;
                }
                std::cerr << "  Added angles for segment " << segment_name << std::endl;
            }

            // Add dihedrals for this molecule instance
            if (molecule_dihedrals_temp.find(mol_type) != molecule_dihedrals_temp.end()) {
                if (!parse_dihedrals_section(molecule_dihedrals_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_dihedrals_temp[mol_type].front();
                    std::cerr << "Error parsing dihedrals for molecule " << mol_type 
                             << " at " << line.source_file << ":" << line.line_number << std::endl;
                    return false;
                }
                std::cerr << "  Added dihedrals for segment " << segment_name << std::endl;
            }

            // Add impropers for this molecule instance
            if (molecule_impropers_temp.find(mol_type) != molecule_impropers_temp.end()) {
                if (!parse_impropers_section(molecule_impropers_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_impropers_temp[mol_type].front();
                    std::cerr << "Error parsing impropers for molecule " << mol_type 
                             << " at " << line.source_file << ":" << line.line_number << std::endl;
                    return false;
                }
                std::cerr << "  Added impropers for segment " << segment_name << std::endl;
            }

            // Add CMAPs for this molecule instance
            if (molecule_cmaps_temp.find(mol_type) != molecule_cmaps_temp.end()) {
                if (!parse_cmaps_section(molecule_cmaps_temp[mol_type], topology, atom_offset)) {
                    const auto& line = molecule_cmaps_temp[mol_type].front();
                    std::cerr << "Error parsing CMAPs for molecule " << mol_type 
                             << " at " << line.source_file << ":" << line.line_number << std::endl;
                    return false;
                }
                std::cerr << "  Added CMAPs for segment " << segment_name << std::endl;
            }

            // Update atom offset for next molecule instance
            atom_offset += atoms_per_molecule[mol_type];
        }
    }

    std::cerr << "\n=== Final topology statistics ===" << std::endl;
    std::cerr << "Total atoms: " << topology.get_num_atoms() << std::endl;
    std::cerr << "Total bonds: " << topology.get_num_bonds() << std::endl;
    std::cerr << "Total angles: " << topology.get_num_angles() << std::endl;
    std::cerr << "Total dihedrals: " << topology.get_num_dihedrals() << std::endl;
    std::cerr << "Total impropers: " << topology.get_num_impropers() << std::endl;
    std::cerr << "Total residues: " << topology.get_num_residues() << std::endl;
    std::cerr << "Total segments: " << topology.get_num_segments() << std::endl;

    return true;
}

bool TopParser::collect_all_lines(const std::string& filename, std::vector<LineInfo>& all_lines,
                                PreprocessorState& pp_state, bool is_main_file) {
    // Check if we've already processed this file
    if (processed_files_.find(filename) != processed_files_.end()) {
        return true;
    }
    processed_files_.insert(filename);

    std::ifstream file(filename);
    if (!file.is_open()) {
        if (is_main_file) {
            std::cerr << "Error: Cannot open main topology file '" << filename << "'" << std::endl;
            return false;
        }
        std::cerr << "Warning: Cannot open included file '" << filename << "'" << std::endl;
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
            std::cerr << "Skipping line due to preprocessor: " << line << std::endl;
            continue;
        }

        // Add normal content line with source information
        all_lines.emplace_back(line, filename, line_number);
    }

    return true;
}

bool TopParser::process_preprocessor_line(const std::string& line, const std::string& parent_file,
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
                std::cerr << "Warning: Include file not found: " << include_path 
                          << " (referenced from " << parent_file << ":" << line_number << ")" << std::endl;
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
        
        std::cerr << "Processing " << directive << " " << macro_name 
                  << ": defined=" << is_defined << ", skip=" << pp_state.skip_section 
                  << ", stack_size=" << pp_state.ifdef_stack.size() << std::endl;
    }
    else if (directive == "#else") {
        if (pp_state.ifdef_stack.empty()) {
            std::cerr << "Warning: Unmatched #else at " << parent_file << ":" << line_number << std::endl;
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
            
            std::cerr << "Processing #else: skip=" << pp_state.skip_section 
                      << ", stack_size=" << pp_state.ifdef_stack.size() 
                      << ", outer_skip=" << outer_skip << std::endl;
        } else {
            std::cerr << "Warning: Multiple #else directives at the same nesting level at "
                      << parent_file << ":" << line_number << std::endl;
        }
    }
    else if (directive == "#endif") {
        if (pp_state.ifdef_stack.empty()) {
            std::cerr << "Warning: Unmatched #endif at " << parent_file << ":" << line_number << std::endl;
            return false;
        }
        
        pp_state.ifdef_stack.pop_back();
        pp_state.else_encountered.pop_back();
        pp_state.skip_section = pp_state.should_skip();
        
        std::cerr << "Processing #endif: skip=" << pp_state.skip_section 
                  << ", stack_size=" << pp_state.ifdef_stack.size() << std::endl;
    }
    else if (directive == "#define") {
        if (!pp_state.should_skip()) {
            std::string macro_name;
            iss >> macro_name;
            std::string macro_value;
            std::getline(iss, macro_value);
            trim(macro_value);
            pp_state.defines[macro_name] = macro_value;
            std::cerr << "Defined macro: " << macro_name << " = " << macro_value << std::endl;
        }
    }
    else if (directive == "#undef") {
        if (!pp_state.should_skip()) {
            std::string macro_name;
            iss >> macro_name;
            pp_state.defines.erase(macro_name);
            std::cerr << "Undefined macro: " << macro_name << std::endl;
        }
    }
    
    return true;
}

void TopParser::parse_sections(const std::vector<LineInfo>& all_lines,
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

std::string TopParser::remove_comment(const std::string& line) {
    size_t comment_pos = line.find(';');
    if (comment_pos != std::string::npos) {
        return line.substr(0, comment_pos);
    }
    return line;
}

bool TopParser::parse_moleculetype_section(const std::vector<LineInfo>& lines, [[maybe_unused]] model::Topology& topology) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
        if (tokens.size() < 2) continue;

        try {
            current_molecule_type_ = tokens[0];
            current_molecule_nrexcl_ = std::stoi(tokens[1]);
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error parsing moleculetype line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_atoms_section(const std::vector<LineInfo>& lines, model::Topology& topology) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
        if (tokens.size() < 8) continue;  // Need at least 8 columns

        try {
            // Parse atom data
            std::string atom_type = tokens[1];
            int residue_number = std::stoi(tokens[2]);
            std::string residue_name = tokens[3];
            std::string atom_name = tokens[4];
            double charge = std::stod(tokens[6]);
            double mass = std::stod(tokens[7]);

            // Add atom to topology using current_molecule_type_ as segment
            topology.add_atom(
                atom_name,
                atom_type,
                charge,
                mass,
                residue_name,
                residue_number,
                current_molecule_type_
            );
        } catch (const std::exception& e) {
            std::cerr << "Error parsing atom line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_bonds_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
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
            std::cerr << "Warning: Error parsing bond line: " << line << " at " 
                     << line_info.source_file << ":" << line_info.line_number 
                     << " - " << e.what() << std::endl;
            continue;  // Continue with next line instead of failing
        }
    }
    return true;  // Return true if we processed all lines (even with some warnings)
}

bool TopParser::parse_angles_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
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
            std::cerr << "Warning: Error parsing angle line: " << line << " at "
                     << line_info.source_file << ":" << line_info.line_number 
                     << " - " << e.what() << std::endl;
            continue;  // Continue with next line instead of failing
        }
    }
    return true;
}

bool TopParser::parse_dihedrals_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    std::cerr << "Parsing " << lines.size() << " dihedral lines" << std::endl;
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
            std::cerr << "Error parsing dihedral line: " << line << " - " << e.what() << std::endl;
            return false;
        }
    }
    
    std::cerr << "Added " << proper_count << " proper dihedrals and " 
              << improper_count << " improper dihedrals" << std::endl;
    return true;
}

bool TopParser::parse_impropers_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
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
            std::cerr << "Error parsing improper line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_molecules_section(const std::vector<LineInfo>& lines, [[maybe_unused]] model::Topology& topology) {
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
            std::cerr << "Error parsing molecules line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

std::string TopParser::trim(std::string& str) {
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

std::vector<std::string> TopParser::split(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    
    while (iss >> token) {
        if (token[0] == ';') break;  // Stop at comments
        tokens.push_back(token);
    }
    
    return tokens;
}

std::string TopParser::resolve_include_path(const std::string& include_path, const std::string& parent_file) {
    namespace fs = std::filesystem;
    
    // Convert parent path to absolute and get its directory
    fs::path parent_path = fs::absolute(parent_file);
    fs::path parent_dir = parent_path.parent_path();
    
    // Debug output
    std::cerr << "Resolving include path: " << include_path << "\n"
              << "Parent file: " << parent_file << "\n"
              << "Parent dir: " << parent_dir.string() << std::endl;
    
    // Simply combine parent directory with include path
    fs::path resolved = parent_dir / include_path;
    if (fs::exists(resolved)) {
        std::cerr << "Found include file at: " << resolved.string() << std::endl;
        return resolved.string();
    }
    
    // If not found, provide error message
    std::cerr << "Warning: Include file not found: " << include_path << "\n"
              << "Tried path: " << resolved.string() << std::endl;
    
    return "";
}

bool TopParser::parse_cmaps_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset) {
    for (const auto& line_info : lines) {
        const std::string& line = line_info.content;
        auto tokens = split(remove_comment(line));
        if (tokens.size() < 6) {  // Need 5 atoms + function type
            std::cerr << "Warning: Skipping CMAP line with insufficient tokens: " << line << std::endl;
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
            
            std::cerr << "Added CMAP between atoms: ";
            for (int idx : cmap_atoms) {
                std::cerr << idx << " ";
            }
            std::cerr << " (function type " << function_type << ")" << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error parsing CMAP line: " << line << " at " 
                     << line_info.source_file << ":" << line_info.line_number 
                     << " - " << e.what() << std::endl;
            return false;
        }
    }
    return true;
}

} // namespace io
} // namespace pygcmc





