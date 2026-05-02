// src/io/topology/topParserPreprocessor.hpp

#pragma once

#include "topParserStructures.hpp"
#include <iostream>
#include <set>

namespace pygcmc {
namespace io {

/**
 * @brief Preprocessor handling functions for TOP file parsing
 */
class TopParserPreprocessor {
public:
    // Debug output control
    static bool& getDebugFlag();

    // Debug printing function
    template<typename... Args>
    static void debug_print(Args&&... args) {
        if (getDebugFlag()) {
            (std::cerr << ... << std::forward<Args>(args));
        }
    }

    // Preprocessor and file handling functions
    static bool collect_all_lines(const std::string& filename, std::vector<LineInfo>& all_lines,
                                 PreprocessorState& pp_state, bool is_main_file,
                                 std::set<std::string>& processed_files);

    static bool process_preprocessor_line(const std::string& line, const std::string& parent_file,
                                        std::vector<LineInfo>& all_lines, PreprocessorState& pp_state,
                                        int line_number, std::set<std::string>& processed_files);

    static void parse_sections(const std::vector<LineInfo>& all_lines,
                              std::map<std::string, std::vector<LineInfo>>& sections);

    // Include handling
    static std::string resolve_include_path(const std::string& include_path, const std::string& parent_file);
};

} // namespace io
} // namespace pygcmc
