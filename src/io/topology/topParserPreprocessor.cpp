// src/io/topology/topParserPreprocessor.cpp

#include "topParserPreprocessor.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <filesystem>

namespace pygcmc {
namespace io {

// Static debug flag definition
static bool debug_enabled_ = false;

bool& TopParserPreprocessor::getDebugFlag() {
    return debug_enabled_;
}

bool TopParserPreprocessor::collect_all_lines(const std::string& filename, std::vector<LineInfo>& all_lines,
                                            PreprocessorState& pp_state, bool is_main_file, 
                                            std::set<std::string>& processed_files) {
    // Check if we've already processed this file
    if (processed_files.find(filename) != processed_files.end()) {
        return true;
    }
    processed_files.insert(filename);

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
        // Simple trim implementation
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), line.end());
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == ';') {
            continue;
        }

        // Handle preprocessor directives
        if (line[0] == '#') {
            if (!process_preprocessor_line(line, filename, all_lines, pp_state, line_number, processed_files)) {
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

bool TopParserPreprocessor::process_preprocessor_line(const std::string& line, const std::string& parent_file,
                                                    std::vector<LineInfo>& all_lines, PreprocessorState& pp_state,
                                                    int line_number, std::set<std::string>& processed_files) {
    std::istringstream iss(line);
    std::string directive;
    iss >> directive;

    if (directive == "#include") {
        // Only process include if we're not in a false ifdef block
        if (!pp_state.should_skip()) {
            std::string include_path;
            std::getline(iss, include_path);
            // Simple trim
            include_path.erase(include_path.begin(), std::find_if(include_path.begin(), include_path.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            include_path.erase(std::find_if(include_path.rbegin(), include_path.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), include_path.end());
            
            // Remove quotes if present
            if (include_path.front() == '"' && include_path.back() == '"') {
                include_path = include_path.substr(1, include_path.length() - 2);
            }

            std::string resolved_path = resolve_include_path(include_path, parent_file);
            if (!resolved_path.empty()) {
                if (!collect_all_lines(resolved_path, all_lines, pp_state, false, processed_files)) {
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
            // Simple trim
            macro_value.erase(macro_value.begin(), std::find_if(macro_value.begin(), macro_value.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            macro_value.erase(std::find_if(macro_value.rbegin(), macro_value.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), macro_value.end());
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

void TopParserPreprocessor::parse_sections(const std::vector<LineInfo>& all_lines,
                                         std::map<std::string, std::vector<LineInfo>>& sections) {
    std::string current_section;
    for (const auto& line_info : all_lines) {
        const std::string& line = line_info.content;
        
        // Check for section header
        if (line[0] == '[') {
            std::string tmp = line;
            // Simple trim
            tmp.erase(tmp.begin(), std::find_if(tmp.begin(), tmp.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            tmp.erase(std::find_if(tmp.rbegin(), tmp.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), tmp.end());
            // Remove brackets and trim again
            if (tmp.front() == '[') tmp.erase(tmp.begin());
            if (!tmp.empty() && tmp.back() == ']') tmp.pop_back();
            // Trim again
            tmp.erase(tmp.begin(), std::find_if(tmp.begin(), tmp.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            tmp.erase(std::find_if(tmp.rbegin(), tmp.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), tmp.end());
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

std::string TopParserPreprocessor::resolve_include_path(const std::string& include_path, const std::string& parent_file) {
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

} // namespace io
} // namespace pygcmc