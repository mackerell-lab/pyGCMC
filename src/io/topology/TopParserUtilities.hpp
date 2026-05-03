// src/io/topology/TopParserUtilities.hpp

#pragma once

#include "TopParserStructures.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

/**
 * @brief Basic utility functions for TOP file parsing
 */
class TopParserUtilities {
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

    // String utility functions
    static std::string remove_comment(const std::string& line);
    static std::string trim(const std::string& str);
    static std::vector<std::string> split(const std::string& str);

    // Type-specific helpers
    static double default_mass_for_atom_type(const std::string& atom_type);
};

} // namespace io
} // namespace pygcmc
