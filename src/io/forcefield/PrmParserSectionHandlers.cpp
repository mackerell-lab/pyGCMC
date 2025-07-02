// src/io/forcefield/PrmParserSectionHandlers.cpp

#include "PrmParserSectionHandlers.hpp"
#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

void PrmParserSectionHandlers::parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    std::string line;
    
    if (debug_output) std::cerr << "\n=== Entering ATOMS/MASS section parsing ===" << std::endl;
    
    // Add debug counter
    int lines_processed = 0;
    int atoms_with_alpha = 0;
    
    while (std::getline(input, line)) {
        lines_processed++;
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
        
        // Note: ATOM lines with ALPHA/THOLE are handled in the pre-scan phase
        // parseAtomsSection should only parse MASS entries, not ATOM lines
        
        // Skip other topology lines from STR files
        if (PrmParserStructures::isTopologyLine(fullLine)) {
            if (debug_output) std::cerr << "Skipping topology line in atoms section: " << fullLine << std::endl;
            continue;
        }
        
        // Check for section end
        if (fullLine == "END" || fullLine == "end" || PrmParserStructures::isBondsSection(fullLine) || PrmParserStructures::isAnglesSection(fullLine) || 
            PrmParserStructures::isDihedralsSection(fullLine) || PrmParserStructures::isImproperSection(fullLine) || 
            PrmParserStructures::isNonbondedSection(fullLine) || PrmParserStructures::isNBFixSection(fullLine)) {
            if (debug_output) std::cerr << "Found section end marker: " << fullLine << std::endl;
            // Note: For atoms section, we use the original line length, not fullLine
            // because fullLine may be concatenated from multiple lines
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
    
    // Add debug summary
    if (debug_output) {
        std::cerr << "\n=== Exiting ATOMS/MASS section parsing ===" << std::endl;
        std::cerr << "Lines processed: " << lines_processed << std::endl;
        std::cerr << "ATOM lines with ALPHA/THOLE found: " << atoms_with_alpha << std::endl;
    }
}

void PrmParserSectionHandlers::parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    std::string line;
    
    if (debug_output) std::cerr << "\n=== Entering BONDS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw bond line: [" << line << "]" << std::endl;
        
        // Save original line length before any modifications
        size_t originalLineLength = line.length();
        
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
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in bonds section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isAnglesSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            if (debug_output) std::cerr << "Found section end marker in bonds: " << line << std::endl;
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
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

void PrmParserSectionHandlers::parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    std::string line;
    
    if (debug_output) std::cerr << "\n=== Entering ANGLES section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output && line.find("ATOM") != std::string::npos) {
            std::cerr << "WARNING: ATOM line in angles section: [" << line << "]" << std::endl;
        }
        if (debug_output) std::cerr << "Raw angle line: [" << line << "]" << std::endl;
        
        // Save original line length before any modifications
        size_t originalLineLength = line.length();
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in angles section" << std::endl;
            continue;
        }
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Skip ATOM lines (they're not angle parameters)
        if (line.find("ATOM ") == 0) {
            if (debug_output) std::cerr << "Skipping ATOM line in angles section: " << line << std::endl;
            continue;
        }
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in angles section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
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