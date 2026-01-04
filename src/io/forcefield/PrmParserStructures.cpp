// src/io/forcefield/PrmParserStructures.cpp

#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>

// Forward declaration for debug access
namespace pygcmc { namespace io { class PRMParser; } }

namespace pygcmc {
namespace io {

// Forward declaration for access to PRMParser debug flag
namespace pygcmc { namespace io { class PRMParser; } }

double PrmParserStructures::safe_stod(const std::string& str, const std::string& context) {
    try {
        return std::stod(str);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to convert '" + str + "' to double in " + context);
    }
}

int PrmParserStructures::safe_stoi(const std::string& str, const std::string& context) {
    try {
        return std::stoi(str);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to convert '" + str + "' to integer in " + context);
    }
}

std::vector<std::string> PrmParserStructures::tokenize(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    
    while (iss >> token) {
        if (token[0] == '!' || token[0] == '#') break;  // Stop at comments
        if (token == "-") continue;  // Skip continuation character
        tokens.push_back(token);
    }
    
    return tokens;
}

bool PrmParserStructures::isCommentLine(const std::string& line) {
    return line.empty() || line[0] == '!' || line[0] == '*' || line[0] == '#';
}

std::string PrmParserStructures::removeComments(const std::string& line) {
    size_t commentPos = line.find_first_of("!*#");
    if (commentPos != std::string::npos) {
        return line.substr(0, commentPos);
    }
    return line;
}

std::string PrmParserStructures::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::string PrmParserStructures::readContinuationLine(std::istream& input, std::string firstLine) {
    // Debug output removed from utility functions for better separation
    std::string fullLine = firstLine;
    std::string currentLine;
    
    bool hasContinuation = false;
    if (!fullLine.empty() && fullLine.back() == '-') {
        hasContinuation = true;
        fullLine.pop_back();
        fullLine = trim(fullLine);
    }
    
    while (hasContinuation) {
        if (!std::getline(input, currentLine)) {
            break;
        }
        // Debug output removed
        
        while (isCommentLine(currentLine)) {
            if (!std::getline(input, currentLine)) {
                return fullLine;
            }
            // Debug output removed
        }
        
        currentLine = removeComments(currentLine);
        currentLine = trim(currentLine);
        
        if (currentLine.empty()) {
            break;
        }
        
        hasContinuation = false;
        if (!currentLine.empty() && currentLine.back() == '-') {
            hasContinuation = true;
            currentLine.pop_back();
            currentLine = trim(currentLine);
        }
        
        if (!fullLine.empty() && !currentLine.empty()) {
            fullLine += " ";
        }
        fullLine += currentLine;
        
        // Debug output removed
    }
    
    // Debug output removed
    return fullLine;
}

bool PrmParserStructures::isAtomsSection(const std::string& line) {
    return line.find("ATOMS") != std::string::npos || line.find("MASS") != std::string::npos;
}

bool PrmParserStructures::isBondsSection(const std::string& line) {
    // Only treat a line that is exactly "BONDS" (after trimming)
    // as the start of the parameter BONDS section.
    const std::string trimmed = trim(line);
    return trimmed == "BONDS";
}

bool PrmParserStructures::isAnglesSection(const std::string& line) {
    // Only treat a line that is exactly "ANGLES" (after trimming) 
    // as the start of the parameter ANGLES section.
    const std::string trimmed = trim(line);
    return trimmed == "ANGLES";
}

bool PrmParserStructures::isDihedralsSection(const std::string& line) {
    // Only treat a line that is exactly "DIHEDRALS" (after trimming)
    // as the start of the parameter DIHEDRALS section.
    const std::string trimmed = trim(line);
    return trimmed == "DIHEDRALS";
}

bool PrmParserStructures::isImproperSection(const std::string& line) {
    return line.find("IMPROPER") != std::string::npos;
}

bool PrmParserStructures::isNonbondedSection(const std::string& line) {
    return line.find("NONBONDED") != std::string::npos || 
           line.find("cutnb") != std::string::npos;
}

bool PrmParserStructures::isNBFixSection(const std::string& line) {
    return line.find("NBFIX") != std::string::npos;
}

bool PrmParserStructures::isAlphaTHoleSection(const std::string& line) {
    const std::string trimmed = trim(line);
    if (trimmed.rfind("ALPHA", 0) == 0) {
        // CHARMM parameter files often start the section with a bare "ALPHA" header
        // (THOLE is then given per line), while some variants use "ALPHA THOLE".
        return true;
    }
    return trimmed.find("ALPHA") != std::string::npos && trimmed.find("THOLE") != std::string::npos;
}

bool PrmParserStructures::isLonePairSection(const std::string& line) {
    return line.find("LONEPAIR") != std::string::npos;
}

bool PrmParserStructures::isAnisotropySection(const std::string& line) {
    return line.find("ANISOTROPY") != std::string::npos;
}

bool PrmParserStructures::isTHoleSection(const std::string& line) {
    return line.find("THOLE") != std::string::npos && line.find("TCUT") != std::string::npos;
}

std::pair<std::string, std::string> PrmParserStructures::make_type_pair(
    const std::string& type1, const std::string& type2) {
    // In CHARMM parameter files, bond parameters are stored in the order they appear
    // Do not reorder them based on string comparison
    return std::make_pair(type1, type2);
}

std::tuple<std::string, std::string, std::string> PrmParserStructures::make_type_triple(
    const std::string& type1, const std::string& type2, const std::string& type3) {
    // For angle parameters in CHARMM force field:
    // 1. The middle atom (type2) must stay in the middle
    // 2. Store parameters in the order they appear in the parameter file
    // This ensures we store the parameters exactly as they appear in the force field
    return std::make_tuple(type1, type2, type3);
}

std::tuple<std::string, std::string, std::string, std::string> PrmParserStructures::make_type_quad(
    const std::string& type1, const std::string& type2, 
    const std::string& type3, const std::string& type4) {
    return std::make_tuple(type1, type2, type3, type4);
}

bool PrmParserStructures::isTopologyLine(const std::string& line) {
    // Check if this line is a topology definition from STR files
    // These should be skipped when parsing parameters
    
    // SPECIAL CASE: ATOM lines with ALPHA/THOLE are parameter definitions, not topology
    if (line.find("ATOM ") == 0 && 
        (line.find("ALPHA") != std::string::npos || line.find("THOLE") != std::string::npos)) {
        return false;  // This is a parameter line, not topology
    }
    
    // Residue and patch definitions
    if (line.find("RESI ") == 0 || line.find("PRES ") == 0) {
        return true;
    }
    
    // Atom definitions within topology (without ALPHA/THOLE)
    if (line.find("ATOM ") == 0) {
        return true;
    }
    
    // Group definitions
    if (line.find("GROUP") == 0) {
        return true;
    }
    
    // Topology bonds (note the space after BOND to distinguish from BONDS section)
    if (line.find("BOND ") == 0) {
        return true;
    }
    
    // Topology impropers (note: IMPR with space, not IMPROPER section)
    if (line.find("IMPR ") == 0) {
        return true;
    }
    
    // Topology dihedrals
    if (line.find("DIHE ") == 0) {
        return true;
    }
    
    // Other topology-specific keywords
    if (line.find("DONOR ") == 0 || line.find("ACCEPTOR ") == 0) {
        return true;
    }
    
    // IC (internal coordinate) definitions
    if (line.find("IC ") == 0) {
        return true;
    }
    
    // PATCH applications
    if (line.find("PATCH") == 0) {
        return true;
    }
    
    // patch first none last none (lowercase patch)
    if (line.find("patch ") == 0) {
        return true;
    }
    
    // Other STR-specific keywords that should be ignored
    if (line.find("NOANG") == 0 || line.find("NODIHE") == 0) {
        return true;
    }
    
    // LONEPAIR definitions (Drude-specific) - but not if they have parameters
    if (line.find("LONEPAIR") == 0) {
        // If line contains parameter keywords, it's a parameter line not topology
        if (line.find("distance") != std::string::npos || 
            line.find("angle") != std::string::npos ||
            line.find("dihe") != std::string::npos) {
            return false;
        }
        return true;
    }
    
    // ANISOTROPY definitions (Drude-specific) - but not if they have parameters
    if (line.find("ANISOTROPY") == 0) {
        // If line contains A11, A22, A33, it's a parameter line not topology
        if (line.find("A11") != std::string::npos || 
            line.find("A22") != std::string::npos ||
            line.find("A33") != std::string::npos) {
            return false;
        }
        return true;
    }
    
    // CMAP definitions in topology
    if (line.find("CMAP") == 0) {
        return true;
    }
    
    // Other topology directives from STR files
    if (line.find("AUTOGENERATE") == 0 ||
        line.find("DECL") == 0 ||
        line.find("DEFA") == 0 ||
        line.find("ANGLE ") == 0) {
        return true;
    }
    
    return false;
}

} // namespace io
} // namespace pygcmc
