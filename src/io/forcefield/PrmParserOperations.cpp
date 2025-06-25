// src/io/forcefield/PrmParserOperations.cpp

#include "PrmParserOperations.hpp"
#include "PrmParserSections.hpp"
#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"
#include <fstream>

namespace pygcmc {
namespace io {

// Static debug flag shared across all parser components
static bool debug_output = false;

bool& PrmParserOperations::getDebugFlag() {
    return debug_output;
}

void PrmParserOperations::parseStream(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    bool inSection = false;
    std::string currentSection;
    
    // Sync debug flags
    PrmParserSections::getDebugFlag() = debug_output;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line" << std::endl;
            continue;
        }
        
        std::string cleanLine = PrmParserStructures::removeComments(line);
        cleanLine = PrmParserStructures::trim(cleanLine);
        if (debug_output) std::cerr << "Cleaned line: [" << cleanLine << "]" << std::endl;
        
        if (cleanLine.empty()) continue;
        
        if (cleanLine == "END") {
            inSection = false;
            currentSection.clear();
            continue;
        }
        
        if (PrmParserStructures::isAtomsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ATOMS section" << std::endl;
            inSection = true;
            currentSection = "ATOMS";
            parseAtomsSection(input, ff);
        } else if (PrmParserStructures::isBondsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found BONDS section" << std::endl;
            inSection = true;
            currentSection = "BONDS";
            parseBondsSection(input, ff);
        } else if (PrmParserStructures::isAnglesSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ANGLES section" << std::endl;
            inSection = true;
            currentSection = "ANGLES";
            parseAnglesSection(input, ff);
        } else if (PrmParserStructures::isDihedralsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found DIHEDRALS section" << std::endl;
            inSection = true;
            currentSection = "DIHEDRALS";
            PrmParserSections::parseDihedralsSection(input, ff);
        } else if (PrmParserStructures::isImproperSection(cleanLine)) {
            if (debug_output) std::cerr << "Found IMPROPER section" << std::endl;
            inSection = true;
            currentSection = "IMPROPER";
            PrmParserSections::parseImproperSection(input, ff);
        } else if (PrmParserStructures::isNonbondedSection(cleanLine)) {
            if (debug_output) std::cerr << "Found NONBONDED section" << std::endl;
            inSection = true;
            currentSection = "NONBONDED";
            
            // If the line starts with cutnb, we need to look for the previous NONBONDED line
            if (cleanLine.find("NONBONDED") == std::string::npos && cleanLine.find("cutnb") != std::string::npos) {
                // Store current position
                auto currentPos = input.tellg();
                std::string prevLine;
                
                // Go back to beginning of file
                input.seekg(0);
                
                // Read lines until we reach our current position
                while (input.tellg() < currentPos && std::getline(input, prevLine)) {
                    if (PrmParserStructures::isCommentLine(prevLine)) continue;
                    prevLine = PrmParserStructures::removeComments(prevLine);
                    prevLine = PrmParserStructures::trim(prevLine);
                    if (prevLine.find("NONBONDED") != std::string::npos) {
                        // Found the NONBONDED line, combine it with current line
                        cleanLine = prevLine + " " + cleanLine;
                        break;
                    }
                }
                
                // Restore position
                input.seekg(currentPos);
            }
            
            PrmParserSections::parseNonbondedSection(input, ff, cleanLine);
        } else if (PrmParserStructures::isNBFixSection(cleanLine)) {
            if (debug_output) std::cerr << "Found NBFIX section" << std::endl;
            inSection = true;
            currentSection = "NBFIX";
            PrmParserSections::parseNBFixSection(input, ff);
        } else if (!inSection) {
            // Handle non-section content if needed
            auto tokens = PrmParserStructures::tokenize(cleanLine);
            if (!tokens.empty()) {
                // Process any standalone parameters or commands
                if (tokens[0] == "set" || tokens[0] == "if" || tokens[0] == "read" || 
                    tokens[0] == "return" || tokens[0] == "BOMLEV" || tokens[0] == "WRNLEV") {
                    // Skip CHARMM control statements
                    continue;
                }
            }
        }
    }
}

void PrmParserOperations::parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering ATOMS/MASS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw atom line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in atoms section" << std::endl;
            continue;
        }
        
        std::string fullLine = PrmParserStructures::readContinuationLine(input, line);
        fullLine = PrmParserStructures::removeComments(fullLine);
        fullLine = PrmParserStructures::trim(fullLine);
        
        if (fullLine.empty()) {
            if (debug_output) std::cerr << "Skipping empty line in atoms section" << std::endl;
            continue;
        }
        
        if (debug_output) std::cerr << "Processed line: [" << fullLine << "]" << std::endl;
        
        // Check for section end
        if (fullLine == "END" || PrmParserStructures::isBondsSection(fullLine) || PrmParserStructures::isAnglesSection(fullLine) || 
            PrmParserStructures::isDihedralsSection(fullLine) || PrmParserStructures::isImproperSection(fullLine) || 
            PrmParserStructures::isNonbondedSection(fullLine) || PrmParserStructures::isNBFixSection(fullLine)) {
            if (debug_output) std::cerr << "Found section end marker: " << fullLine << std::endl;
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(fullLine);
        if (debug_output) {
            std::cerr << "Tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.size() >= 4 && tokens[0] == "MASS") {
            std::string atomType = tokens[2];
            double mass = PrmParserStructures::safe_stod(tokens[3], "atom mass for type " + atomType);
            ff.add_atom_mass(atomType, mass);
            if (debug_output) std::cerr << "Added atom mass: " << atomType << " = " << mass << std::endl;
        } else {
            if (debug_output) std::cerr << "Skipping line: not a valid MASS entry" << std::endl;
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished ATOMS/MASS section parsing ===" << std::endl;
        std::cerr << "Final atom_masses size: " << ff.get_num_atom_types() << std::endl;
        std::cerr << "\nStored atom masses:" << std::endl;
        for (const auto& pair : ff.get_atom_masses()) {
            std::cerr << pair.first << " = " << pair.second << std::endl;
        }
    }
}

void PrmParserOperations::parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering BONDS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw bond line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in bonds section" << std::endl;
            continue;
        }
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) {
            if (debug_output) std::cerr << "Skipping empty line in bonds section" << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isAnglesSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            if (debug_output) std::cerr << "Found section end marker in bonds: " << line << std::endl;
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (debug_output) {
            std::cerr << "Bond tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.size() >= 4) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            double kb = PrmParserStructures::safe_stod(tokens[2], "bond Kb for " + type1 + "-" + type2);
            double b0 = PrmParserStructures::safe_stod(tokens[3], "bond b0 for " + type1 + "-" + type2);
            
            if (debug_output) {
                std::cerr << "\nProcessing bond: " << type1 << "-" << type2 << std::endl;
                std::cerr << "  kb = " << kb << ", b0 = " << b0 << std::endl;
            }
            
            ff.add_bond_params(type1, type2, kb, b0);
            
            if (debug_output) {
                std::cerr << "  Current bond_params size: " << ff.get_num_bond_types() << std::endl;
            }
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished BONDS section parsing ===" << std::endl;
        std::cerr << "Final bond_params size: " << ff.get_num_bond_types() << std::endl;
    }
}

void PrmParserOperations::parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering ANGLES section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw angle line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in angles section" << std::endl;
            continue;
        }
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (tokens.size() >= 5) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            std::string type3 = tokens[2];
            double ktheta = PrmParserStructures::safe_stod(tokens[3], "angle Ktheta for " + type1 + "-" + type2 + "-" + type3);
            double theta0 = PrmParserStructures::safe_stod(tokens[4], "angle theta0 for " + type1 + "-" + type2 + "-" + type3);
            
            double kub = 0.0;
            double s0 = 0.0;
            if (tokens.size() >= 7) {
                kub = PrmParserStructures::safe_stod(tokens[5], "angle Kub for " + type1 + "-" + type2 + "-" + type3);
                s0 = PrmParserStructures::safe_stod(tokens[6], "angle S0 for " + type1 + "-" + type2 + "-" + type3);
            }
            
            if (debug_output) {
                std::cerr << "\nProcessing angle: " << type1 << "-" << type2 << "-" << type3 << std::endl;
                std::cerr << "  ktheta = " << ktheta << ", theta0 = " << theta0 << std::endl;
            }
            
            ff.add_angle_params(type1, type2, type3, ktheta, theta0, kub, s0);
            
            if (debug_output) {
                std::cerr << "  Current angle_params size: " << ff.get_num_angle_types() << std::endl;
            }
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished ANGLES section parsing ===" << std::endl;
        std::cerr << "Final angle_params size: " << ff.get_num_angle_types() << std::endl;
    }
}

} // namespace io
} // namespace pygcmc