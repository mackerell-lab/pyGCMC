// src/io/topology/TopParserUtilities.cpp

#include "TopParserUtilities.hpp"
#include <sstream>
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace io {

// Static debug flag definition
static bool debug_enabled_ = false;

bool& TopParserUtilities::getDebugFlag() {
    return debug_enabled_;
}

std::string TopParserUtilities::remove_comment(const std::string& line) {
    size_t comment_pos = line.find(';');
    if (comment_pos != std::string::npos) {
        return line.substr(0, comment_pos);
    }
    return line;
}

std::string TopParserUtilities::trim(const std::string& str) {
    std::string trimmed = str;

    // Trim leading spaces
    trimmed.erase(trimmed.begin(), std::find_if(trimmed.begin(), trimmed.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));

    // Trim trailing spaces
    trimmed.erase(std::find_if(trimmed.rbegin(), trimmed.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), trimmed.end());

    return trimmed;
}

std::vector<std::string> TopParserUtilities::split(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;

    while (iss >> token) {
        if (token[0] == ';') break;  // Stop at comments
        tokens.push_back(token);
    }

    return tokens;
}

double TopParserUtilities::default_mass_for_atom_type(const std::string& atom_type) {
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
