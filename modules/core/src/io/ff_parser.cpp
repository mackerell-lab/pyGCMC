// modules/core/src/io/ff_parser.cpp

#include "pygcmc/core/io/ff_parser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <algorithm>

namespace pygcmc {
namespace core {
namespace io {

bool FFParser::parse(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open force field file: " << filename << std::endl;
        return false;
    }

    // Clear existing data
    nonbonded_params_.clear();
    nbfix_params_.clear();

    std::string line;
    while (std::getline(ifs, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Remove leading/trailing whitespace
        auto startPos = line.find_first_not_of(" \t\r\n");
        if (startPos == std::string::npos) continue;
        line.erase(0, startPos);
        auto endPos = line.find_last_not_of(" \t\r\n");
        if (endPos != std::string::npos) {
            line.erase(endPos + 1);
        }

        // Convert to uppercase for keyword matching
        std::string uline = line;
        std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);

        // Parse sections based on keywords
        if (uline.rfind("NONBONDED", 0) == 0) {
            parse_nonbonded_section(ifs);
        }
        else if (uline.rfind("NBFIX", 0) == 0) {
            parse_nbfix_section(ifs);
        }
        else if (uline == "END") {
            break;
        }
    }

    ifs.close();
    return true;
}

void FFParser::parse_nonbonded_section(std::istream& in) {
    std::string line;
    // Skip header lines until we find actual data
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        
        // Remove comments
        auto cpos = line.find('!');
        if (cpos != std::string::npos) {
            line = line.substr(0, cpos);
        }

        // Remove leading/trailing whitespace
        auto startPos = line.find_first_not_of(" \t\r\n");
        if (startPos == std::string::npos) continue;
        line.erase(0, startPos);
        auto endPos = line.find_last_not_of(" \t\r\n");
        if (endPos != std::string::npos) {
            line.erase(endPos + 1);
        }

        // Check for section end or next section
        {
            std::string tmp = line;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            if (tmp.rfind("NBFIX", 0) == 0 || tmp == "END") {
                in.seekg(-static_cast<int>(line.size())-1, std::ios::cur);
                return;
            }
        }

        // Parse nonbonded parameters
        // Format: atomType ignored epsilon Rmin/2 [ignored ignored ignored]
        std::istringstream iss(line);
        std::string atomType;
        double ignored, eps, rmin;
        
        if (!(iss >> atomType)) continue;
        if (atomType[0] == '!') continue;  // Skip comment lines
        
        if (!(iss >> ignored >> eps >> rmin)) {
            continue;
        }

        // Store parameters (use absolute value of epsilon)
        nonbonded_params_[atomType] = ForceFieldPair(std::abs(eps), rmin);
    }
}

void FFParser::parse_nbfix_section(std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // Remove comments
        auto cpos = line.find('!');
        if (cpos != std::string::npos) {
            line = line.substr(0, cpos);
        }

        // Remove leading/trailing whitespace
        auto startPos = line.find_first_not_of(" \t\r\n");
        if (startPos == std::string::npos) continue;
        line.erase(0, startPos);
        auto endPos = line.find_last_not_of(" \t\r\n");
        if (endPos != std::string::npos) {
            line.erase(endPos + 1);
        }

        // Check for section end or next section
        {
            std::string tmp = line;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            if (tmp == "END" || tmp.rfind("NONBONDED", 0) == 0) {
                in.seekg(-static_cast<int>(line.size())-1, std::ios::cur);
                return;
            }
        }

        // Parse NBFIX parameters
        // Format: type1 type2 epsilon Rmin
        std::istringstream iss(line);
        std::string t1, t2;
        double eps, rmin;
        
        if (!(iss >> t1)) continue;
        if (t1[0] == '!') continue;  // Skip comment lines
        
        if (!(iss >> t2 >> eps >> rmin)) {
            continue;
        }

        // Store parameters (both directions due to symmetry)
        // Use absolute value of epsilon
        ForceFieldPair ff(std::abs(eps), rmin);
        nbfix_params_[std::make_pair(t1, t2)] = ff;
        nbfix_params_[std::make_pair(t2, t1)] = ff;
    }
}

} // namespace io
} // namespace core
} // namespace pygcmc

