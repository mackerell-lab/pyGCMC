// src/io/topology/TopParserStructures.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <string>
#include <vector>
#include <map>

namespace pygcmc {
namespace io {

/**
 * @brief Structure to track line source information
 */
struct LineInfo {
    std::string content;      ///< The actual line content
    std::string source_file;  ///< Source file path
    int line_number;          ///< Line number in source file

    LineInfo(const std::string& content, const std::string& file, int line)
        : content(content), source_file(file), line_number(line) {}
};

/**
 * @brief Preprocessor state for handling #ifdef, #define, etc.
 */
struct PreprocessorState {
    std::map<std::string, std::string> defines;  ///< #define macros
    std::vector<bool> ifdef_stack;               ///< Stack for #ifdef/#ifndef nesting
    std::vector<bool> else_encountered;          ///< Track if #else was encountered at each nesting level
    bool skip_section = false;                   ///< Whether to skip current section due to #ifdef

    bool should_skip() const {
        // Skip if any level in the stack is false
        for (bool val : ifdef_stack) {
            if (!val) return true;
        }
        return false;
    }
};

} // namespace io
} // namespace pygcmc
