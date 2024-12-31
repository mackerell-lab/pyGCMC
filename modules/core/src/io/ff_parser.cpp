// modules/core/src/io/ff_parser.cpp

#include "pygcmc/core/io/ff_parser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <cmath>

namespace pygcmc {
namespace core {
namespace io {

bool FFParser::parse(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open force field file: " << filename << std::endl;
        return false;
    }

    std::string line;
    bool inNonbondedSection = false;
    bool inNBFIXSection = false;

    while (std::getline(ifs, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Remove comments and trim whitespace
        auto commentPos = line.find('!');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        auto startPos = line.find_first_not_of(" \t\r\n");
        if (startPos == std::string::npos) continue;
        line.erase(0, startPos);
        auto endPos = line.find_last_not_of(" \t\r\n");
        if (endPos != std::string::npos) {
            line.erase(endPos + 1);
        }

        // Skip empty lines after removing comments
        if (line.empty()) continue;

        // Convert to uppercase for keyword comparison
        std::string uline = line;
        std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);

        // Detect sections
        if (uline.rfind("NONBONDED", 0) == 0) {
            inNonbondedSection = true;
            inNBFIXSection = false;
            continue;
        }
        else if (uline.rfind("NBFIX", 0) == 0) {
            inNonbondedSection = false;
            inNBFIXSection = true;
            continue;
        }
        else if (uline == "END") {
            inNonbondedSection = false;
            inNBFIXSection = false;
            continue;
        }

        // Parse lines based on current section
        if (inNonbondedSection) {
            parse_nonbonded_line(line);
        }
        else if (inNBFIXSection) {
            parse_nbfix_line(line);
        }
    }

    ifs.close();
    return true;
}

void FFParser::parse_nonbonded_line(const std::string& line) {
    // Skip header lines
    if (line.find("nbxmod") != std::string::npos ||
        line.find("cutnb") != std::string::npos ||
        line.find("ctofnb") != std::string::npos ||
        line.find("ctonnb") != std::string::npos ||
        line.find("eps") != std::string::npos ||
        line.find("e14fac") != std::string::npos ||
        line.find("wmin") != std::string::npos ||
        line.find("NONBONDED") != std::string::npos) {
        return;
    }

    // Skip lines that start with special characters
    if (line[0] == '-' || line[0] == '*' || line[0] == '@' || line[0] == '#' || line[0] == '!') {
        return;
    }

    // Parse atom parameters
    std::istringstream iss(line);
    std::string atomType;
    double ignore1, epsilon, rminHalf;
    double eps14 = 0.0, rmin14Half = 0.0;
    
    if (iss >> atomType >> ignore1 >> epsilon >> rminHalf) {
        // Try to read 1-4 parameters if present
        double ignore2;
        if (iss >> ignore2 >> eps14 >> rmin14Half) {
            // Store 1-4 parameters if needed
        }

        // Process standard parameters:
        // epsilon is negative in file (e.g. -0.2), take absolute value
        // rminHalf is Rmin/2, multiply by 2 to get Rmin
        double rmin = rminHalf * 2.0;
        nonbonded_params_[atomType] = ForceFieldPair(rmin, std::fabs(epsilon));

        // Debug output
        // std::cout << "Parsed nonbonded: " << atomType 
        //           << " | epsilon=" << std::fabs(epsilon) 
        //           << ", rmin=" << rmin << std::endl;
    }
}

void FFParser::parse_nbfix_line(const std::string& line) {
    // Skip header lines
    if (line.find("Emin") != std::string::npos || 
        line.find("kcal/mol") != std::string::npos ||
        line.find("NBFIX") != std::string::npos) {
        return;
    }

    // Skip lines that start with special characters
    if (line[0] == '-' || line[0] == '*' || line[0] == '@' || line[0] == '#' || line[0] == '!') {
        return;
    }

    // Parse NBFIX line
    std::istringstream iss(line);
    std::string atom_type1, atom_type2;
    double epsilon, rmin;

    if (iss >> atom_type1 >> atom_type2 >> epsilon >> rmin) {
        // In CHARMM format, epsilon is negative
        auto key = std::make_pair(atom_type1, atom_type2);
        nbfix_params_[key] = ForceFieldPair(rmin, std::abs(epsilon));
        
        // Add reverse pair
        auto key_rev = std::make_pair(atom_type2, atom_type1);
        nbfix_params_[key_rev] = ForceFieldPair(rmin, std::abs(epsilon));
        
        // Debug output
        // std::cout << "Parsed NBFIX: " << atom_type1 << " - " << atom_type2 
        //           << " | epsilon=" << std::abs(epsilon) 
        //           << ", rmin=" << rmin << std::endl;
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
            // 仅在参数未设置（为 nan）时更新
            if (std::isnan(atom.forcefield_epsilon) || std::isnan(atom.forcefield_rmin)) {
                atom.forcefield_epsilon = it->second.epsilon;
                atom.forcefield_rmin = it->second.rmin;
                updated++;
            }
        }
    }
    return updated;
}

int FFParser::update_pdb_atoms(std::vector<PDBAtom*>& atoms) const {
    int updated = 0;
    for (auto* atom : atoms) {
        if (atom->topo_type.empty()) {
            continue;
        }
        auto it = nonbonded_params_.find(atom->topo_type);
        if (it != nonbonded_params_.end()) {
            // 仅在参数未设置（为 nan）时更新
            if (std::isnan(atom->forcefield_epsilon) || std::isnan(atom->forcefield_rmin)) {
                atom->forcefield_epsilon = it->second.epsilon;
                atom->forcefield_rmin = it->second.rmin;
                updated++;
            }
        }
    }
    return updated;
}

} // namespace io
} // namespace core
} // namespace pygcmc
