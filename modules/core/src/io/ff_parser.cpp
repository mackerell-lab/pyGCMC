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
    bool inHeader = true;

    while (std::getline(in, line)) {
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

        // Check for section endings (convert to uppercase for comparison)
        std::string uline = line;
        std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);
        if (uline.rfind("NBFIX", 0) == 0) {
            parse_nbfix_section(in);
            return;
        }
        if (uline == "END") {
            return;
        }

        // Skip comments and continuation lines
        if (line[0] == '!' || line[0] == '-') continue;

        // Skip header lines until we find actual data
        if (inHeader) {
            // Check if this line contains the nbxmod keyword
            if (line.find("nbxmod") != std::string::npos) {
                continue;
            }
            // If we get here and it's not a comment, we're past the header
            if (line[0] != '!') {
                inHeader = false;
            }
            continue;
        }

        // Parse the data line
        std::istringstream iss(line);
        std::string atomType;
        double ignore1, epsilon, rminHalf;
        double eps14 = 0.0, rmin14Half = 0.0;
        
        // Try to read both standard and 1-4 parameters
        // Format: atomType ignored -epsilon Rmin/2 [ignored eps14 Rmin14/2]
        if (iss >> atomType >> ignore1 >> epsilon >> rminHalf) {
            // Skip if this is a comment line that happens to have numbers
            if (atomType[0] == '!') continue;
            
            // Try to read 1-4 parameters if they exist
            double ignore2;
            if (iss >> ignore2 >> eps14 >> rmin14Half) {
                // Store 1-4 parameters if needed
                // nonbonded_params_14_[atomType] = ForceFieldPair(rmin14Half * 2.0, std::fabs(eps14));
            }

            // Process standard parameters:
            // epsilon is negative in file (-0.2), take absolute value (0.2)
            // rminHalf is 1.85, multiply by 2 to get Rmin (3.70)
            double rmin = rminHalf * 2.0;
            nonbonded_params_[atomType] = ForceFieldPair(rmin, std::fabs(epsilon));
        }
    }
}

void FFParser::parse_nbfix_section(std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '!') continue;

        // Check for section endings (convert to uppercase for comparison)
        std::string uline = line;
        std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);
        if (uline.rfind("NONBONDED", 0) == 0) {
            parse_nonbonded_section(in);
            return;
        }
        if (uline == "END") {
            return;
        }

        // Parse the line
        std::istringstream iss(line);
        std::string atom_type1, atom_type2;
        double epsilon, rmin;

        // CHARMM format has: type1 type2 -epsilon Rmin [ignored] [ignored] -epsilon_14 Rmin_14
        // We only care about the first epsilon and Rmin
        if (iss >> atom_type1 >> atom_type2 >> epsilon >> rmin) {
            // In CHARMM format, the values are epsilon and Rmin
            // Note: epsilon is positive in the parameter file
            auto key = std::make_pair(atom_type1, atom_type2);
            nbfix_params_[key] = ForceFieldPair(rmin, std::abs(epsilon));
            
            // Add reverse pair with same parameters
            auto key_rev = std::make_pair(atom_type2, atom_type1);
            nbfix_params_[key_rev] = ForceFieldPair(rmin, std::abs(epsilon));
        }
    }
}

int FFParser::update_pdb_atoms(std::vector<PDBAtom>& atoms) const {
    int updated = 0;
    for (auto& atom : atoms) {
        if (atom.topo_type.empty()) {
            continue;
        }
        auto it = nonbonded_params_.find(atom.topo_type);
        if (it != nonbonded_params_.end()) {
            atom.forcefield_epsilon = it->second.epsilon;
            atom.forcefield_rmin = it->second.rmin;
            updated++;
        }
    }
    return updated;
}

} // namespace io
} // namespace core
} // namespace pygcmc

