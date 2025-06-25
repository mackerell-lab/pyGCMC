// src/io/forcefield/PrmParserOperations.cpp

#include "PrmParserOperations.hpp"
#include "PrmParserMain.hpp"
#include "model/ModelModule.hpp"
#include <fstream>

namespace pygcmc {
namespace io {

void PrmParserOperations::parse_string(const std::string& content, pygcmc::model::ForceField& ff) {
    std::istringstream iss(content);
    parseStream(iss, ff);
}

void PrmParserOperations::parse_file_to_forcefield(const std::string& filename, pygcmc::model::ForceField& ff) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open parameter file: " + filename);
    }
    parseStream(file, ff);
}

pygcmc::model::ForceField PrmParserOperations::parse_file(const std::string& filename) {
    pygcmc::model::ForceField ff;
    parse_file_to_forcefield(filename, ff);
    return ff;
}

pygcmc::model::ForceField PrmParserOperations::parse_files(const std::vector<std::string>& filenames) {
    pygcmc::model::ForceField ff;
    for (const auto& filename : filenames) {
        parse_file_to_forcefield(filename, ff);
    }
    return ff;
}

void PrmParserOperations::parseStream(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    bool inSection = false;
    std::string currentSection;
    
    while (std::getline(input, line)) {
        if (PRMParser::debug_output) std::cerr << "Raw line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (PRMParser::debug_output) std::cerr << "Skipping comment line" << std::endl;
            continue;
        }
        
        std::string cleanLine = PrmParserStructures::removeComments(line);
        cleanLine = PrmParserStructures::trim(cleanLine);
        if (PRMParser::debug_output) std::cerr << "Cleaned line: [" << cleanLine << "]" << std::endl;
        
        if (cleanLine.empty()) continue;
        
        if (cleanLine == "END") {
            inSection = false;
            currentSection.clear();
            continue;
        }
        
        if (PrmParserStructures::isAtomsSection(cleanLine)) {
            if (PRMParser::debug_output) std::cerr << "Found ATOMS section" << std::endl;
            inSection = true;
            currentSection = "ATOMS";
            parseAtomsSection(input, ff);
        } else if (PrmParserStructures::isBondsSection(cleanLine)) {
            if (PRMParser::debug_output) std::cerr << "Found BONDS section" << std::endl;
            inSection = true;
            currentSection = "BONDS";
            parseBondsSection(input, ff);
        } else if (PrmParserStructures::isAnglesSection(cleanLine)) {
            if (PRMParser::debug_output) std::cerr << "Found ANGLES section" << std::endl;
            inSection = true;
            currentSection = "ANGLES";
            parseAnglesSection(input, ff);
        } else if (PrmParserStructures::isDihedralsSection(cleanLine)) {
            if (PRMParser::debug_output) std::cerr << "Found DIHEDRALS section" << std::endl;
            inSection = true;
            currentSection = "DIHEDRALS";
            parseDihedralsSection(input, ff);
        } else if (PrmParserStructures::isImproperSection(cleanLine)) {
            if (PRMParser::debug_output) std::cerr << "Found IMPROPER section" << std::endl;
            inSection = true;
            currentSection = "IMPROPER";
            parseImproperSection(input, ff);
        } else if (PrmParserStructures::isNonbondedSection(cleanLine)) {
            if (PRMParser::debug_output) std::cerr << "Found NONBONDED section" << std::endl;
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
            
            parseNonbondedSection(input, ff, cleanLine);
        } else if (PrmParserStructures::isNBFixSection(cleanLine)) {
            if (PRMParser::debug_output) std::cerr << "Found NBFIX section" << std::endl;
            inSection = true;
            currentSection = "NBFIX";
            parseNBFixSection(input, ff);
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
    if (PRMParser::debug_output) std::cerr << "\n=== Entering ATOMS/MASS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (PRMParser::debug_output) std::cerr << "Raw atom line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (PRMParser::debug_output) std::cerr << "Skipping comment line in atoms section" << std::endl;
            continue;
        }
        
        std::string fullLine = PrmParserStructures::readContinuationLine(input, line);
        fullLine = PrmParserStructures::removeComments(fullLine);
        fullLine = PrmParserStructures::trim(fullLine);
        
        if (fullLine.empty()) {
            if (PRMParser::debug_output) std::cerr << "Skipping empty line in atoms section" << std::endl;
            continue;
        }
        
        if (PRMParser::debug_output) std::cerr << "Processed line: [" << fullLine << "]" << std::endl;
        
        // Check for section end
        if (fullLine == "END" || PrmParserStructures::isBondsSection(fullLine) || PrmParserStructures::isAnglesSection(fullLine) || 
            PrmParserStructures::isDihedralsSection(fullLine) || PrmParserStructures::isImproperSection(fullLine) || 
            PrmParserStructures::isNonbondedSection(fullLine) || PrmParserStructures::isNBFixSection(fullLine)) {
            if (PRMParser::debug_output) std::cerr << "Found section end marker: " << fullLine << std::endl;
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(fullLine);
        if (PRMParser::debug_output) {
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
            if (PRMParser::debug_output) std::cerr << "Added atom mass: " << atomType << " = " << mass << std::endl;
        } else {
            if (PRMParser::debug_output) std::cerr << "Skipping line: not a valid MASS entry" << std::endl;
        }
    }
    
    if (PRMParser::debug_output) {
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
    if (PRMParser::debug_output) std::cerr << "\n=== Entering BONDS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (PRMParser::debug_output) std::cerr << "Raw bond line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (PRMParser::debug_output) std::cerr << "Skipping comment line in bonds section" << std::endl;
            continue;
        }
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) {
            if (PRMParser::debug_output) std::cerr << "Skipping empty line in bonds section" << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isAnglesSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            if (PRMParser::debug_output) std::cerr << "Found section end marker in bonds: " << line << std::endl;
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (PRMParser::debug_output) {
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
            
            if (PRMParser::debug_output) {
                std::cerr << "\nProcessing bond: " << type1 << "-" << type2 << std::endl;
                std::cerr << "  kb = " << kb << ", b0 = " << b0 << std::endl;
            }
            
            ff.add_bond_params(type1, type2, kb, b0);
            
            if (PRMParser::debug_output) {
                std::cerr << "  Current bond_params size: " << ff.get_num_bond_types() << std::endl;
            }
        }
    }
    
    if (PRMParser::debug_output) {
        std::cerr << "\n=== Finished BONDS section parsing ===" << std::endl;
        std::cerr << "Final bond_params size: " << ff.get_num_bond_types() << std::endl;
    }
}

void PrmParserOperations::parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    if (PRMParser::debug_output) std::cerr << "\n=== Entering ANGLES section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (PRMParser::debug_output) std::cerr << "Raw angle line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (PRMParser::debug_output) std::cerr << "Skipping comment line in angles section" << std::endl;
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
            
            if (PRMParser::debug_output) {
                std::cerr << "\nProcessing angle: " << type1 << "-" << type2 << "-" << type3 << std::endl;
                std::cerr << "  ktheta = " << ktheta << ", theta0 = " << theta0 << std::endl;
            }
            
            ff.add_angle_params(type1, type2, type3, ktheta, theta0, kub, s0);
            
            if (PRMParser::debug_output) {
                std::cerr << "  Current angle_params size: " << ff.get_num_angle_types() << std::endl;
            }
        }
    }
    
    if (PRMParser::debug_output) {
        std::cerr << "\n=== Finished ANGLES section parsing ===" << std::endl;
        std::cerr << "Final angle_params size: " << ff.get_num_angle_types() << std::endl;
    }
}

void PrmParserOperations::parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (tokens.size() >= 7) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            std::string type3 = tokens[2];
            std::string type4 = tokens[3];
            double kchi = PrmParserStructures::safe_stod(tokens[4], "dihedral Kchi for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            int n = PrmParserStructures::safe_stoi(tokens[5], "dihedral n for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            double delta = PrmParserStructures::safe_stod(tokens[6], "dihedral delta for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            
            ff.add_dihedral_params(type1, type2, type3, type4, kchi, n, delta);
        }
    }
}

void PrmParserOperations::parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isDihedralsSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (tokens.size() >= 6) {
            try {
                std::string type1 = tokens[0];
                std::string type2 = tokens[1];
                std::string type3 = tokens[2];
                std::string type4 = tokens[3];
                double kpsi = PrmParserStructures::safe_stod(tokens[4], "improper Kpsi");
                double psi0 = PrmParserStructures::safe_stod(tokens[5], "improper psi0");
                
                ff.add_improper_params(type1, type2, type3, type4, kpsi, psi0);
            } catch (const std::exception& e) {
                if (PRMParser::debug_output) std::cerr << "Warning: Skipping improper line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

// Split the large parseNonbondedSection and parseNBFixSection into separate file for better readability
// This implementation keeps them here for completeness but could be further modularized if needed

void PrmParserOperations::parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine) {
    std::string line = firstLine;
    std::string fullLine = PrmParserStructures::readContinuationLine(input, line);
    
    if (PRMParser::debug_output) std::cerr << "=== Entering NONBONDED section parsing ===" << std::endl;
    if (PRMParser::debug_output) std::cerr << "First line: [" << firstLine << "]" << std::endl;
    
    // Parse header parameters
    auto tokens = PrmParserStructures::tokenize(fullLine);
    if (!tokens.empty()) {
        if (PRMParser::debug_output) {
            std::cerr << "Initial tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        // Skip the NONBONDED keyword
        if (tokens[0] == "NONBONDED") {
            if (PRMParser::debug_output) std::cerr << "Skipping NONBONDED keyword" << std::endl;
            tokens.erase(tokens.begin());
        }
        
        // Process parameters
        auto& params = ff.get_nonbonded_params();
        for (size_t i = 0; i < tokens.size(); ++i) {
            if (PRMParser::debug_output) std::cerr << "Processing parameter token: [" << tokens[i] << "]" << std::endl;
            
            if (tokens[i] == "nbxmod" && i + 1 < tokens.size()) {
                params.nbxmod = PrmParserStructures::safe_stoi(tokens[++i], "nbxmod");
                if (PRMParser::debug_output) std::cerr << "Set nbxmod = " << params.nbxmod << std::endl;
            } else if (tokens[i] == "cutnb" && i + 1 < tokens.size()) {
                params.cutnb = PrmParserStructures::safe_stod(tokens[++i], "cutnb");
                if (PRMParser::debug_output) std::cerr << "Set cutnb = " << params.cutnb << std::endl;
            } else if (tokens[i] == "ctofnb" && i + 1 < tokens.size()) {
                params.ctofnb = PrmParserStructures::safe_stod(tokens[++i], "ctofnb");
                if (PRMParser::debug_output) std::cerr << "Set ctofnb = " << params.ctofnb << std::endl;
            } else if (tokens[i] == "ctonnb" && i + 1 < tokens.size()) {
                params.ctonnb = PrmParserStructures::safe_stod(tokens[++i], "ctonnb");
                if (PRMParser::debug_output) std::cerr << "Set ctonnb = " << params.ctonnb << std::endl;
            } else if (tokens[i] == "eps" && i + 1 < tokens.size()) {
                params.eps = PrmParserStructures::safe_stod(tokens[++i], "eps");
                if (PRMParser::debug_output) std::cerr << "Set eps = " << params.eps << std::endl;
            } else if (tokens[i] == "e14fac" && i + 1 < tokens.size()) {
                params.e14fac = PrmParserStructures::safe_stod(tokens[++i], "e14fac");
                if (PRMParser::debug_output) std::cerr << "Set e14fac = " << params.e14fac << std::endl;
            } else if (tokens[i] == "wmin" && i + 1 < tokens.size()) {
                params.wmin = PrmParserStructures::safe_stod(tokens[++i], "wmin");
                if (PRMParser::debug_output) std::cerr << "Set wmin = " << params.wmin << std::endl;
            } else if (tokens[i] == "cdiel") {
                params.cdiel = true;
                if (PRMParser::debug_output) std::cerr << "Set cdiel = true" << std::endl;
            } else if (tokens[i] == "fshift") {
                params.fshift = true;
                if (PRMParser::debug_output) std::cerr << "Set fshift = true" << std::endl;
            } else if (tokens[i] == "vatom") {
                params.vatom = true;
                if (PRMParser::debug_output) std::cerr << "Set vatom = true" << std::endl;
            } else if (tokens[i] == "vdistance") {
                params.vdistance = true;
                if (PRMParser::debug_output) std::cerr << "Set vdistance = true" << std::endl;
            } else if (tokens[i] == "vfswitch") {
                params.vfswitch = true;
                if (PRMParser::debug_output) std::cerr << "Set vfswitch = true" << std::endl;
            }
        }
    }
    
    if (PRMParser::debug_output) std::cerr << "\n=== Starting atom type parameters parsing ===" << std::endl;
    if (PRMParser::debug_output) std::cerr << "Current lj_params map size: " << ff.get_num_lj_params() << std::endl;
    
    // Parse atom type parameters
    while (std::getline(input, line)) {
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) {
            if (PRMParser::debug_output) std::cerr << "Skipping empty line" << std::endl;
            continue;
        }
        
        if (PRMParser::debug_output) std::cerr << "Processing cleaned line: [" << line << "]" << std::endl;
        tokens = PrmParserStructures::tokenize(line);
        if (PRMParser::debug_output) {
            std::cerr << "Tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.empty()) continue;
        
        // Check for section end or new section
        if (line == "END" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isDihedralsSection(line) || 
            PrmParserStructures::isImproperSection(line)) {
            if (PRMParser::debug_output) std::cerr << "Found section end marker: " << tokens[0] << std::endl;
            break;
        }
        
        // If we find NBFIX, parse it as a new section
        if (PrmParserStructures::isNBFixSection(line)) {
            if (PRMParser::debug_output) std::cerr << "Found NBFIX section" << std::endl;
            parseNBFixSection(input, ff);
            break;
        }
        
        // Format: atomType ignored epsilon Rmin [ignored ignored ignored ignored]
        if (tokens.size() >= 4) {
            std::string atomType = tokens[0];
            // Skip the ignored value (usually 0.0)
            try {
                double epsilon = PrmParserStructures::safe_stod(tokens[2], "LJ epsilon for " + atomType);
                double rmin_half = PrmParserStructures::safe_stod(tokens[3], "LJ Rmin/2 for " + atomType);
                
                if (PRMParser::debug_output) {
                    std::cerr << "\n*** Parsing atom type: " << atomType << " ***" << std::endl;
                    std::cerr << "  epsilon = " << epsilon << std::endl;
                    std::cerr << "  rmin_half = " << rmin_half << std::endl;
                }
                
                ff.add_lj_params(atomType, epsilon, rmin_half);
                
                if (PRMParser::debug_output) {
                    std::cerr << "Successfully stored " << atomType << " parameters" << std::endl;
                    std::cerr << "Current lj_params map size: " << ff.get_num_lj_params() << std::endl;
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to parse LJ parameters for " + atomType + ": " + e.what());
            }
        } else if (!tokens.empty() && !PrmParserStructures::isCommentLine(line)) {
            throw std::runtime_error("Malformed NONBONDED parameters in line: " + line + 
                                   "\nExpected at least 4 tokens, got " + std::to_string(tokens.size()));
        }
    }
    
    if (PRMParser::debug_output) std::cerr << "\n=== Finished NONBONDED section parsing ===" << std::endl;
    if (PRMParser::debug_output) std::cerr << "Final lj_params map size: " << ff.get_num_lj_params() << std::endl;
}

void PrmParserOperations::parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Check for section end or new section
        if (line == "END" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isDihedralsSection(line) || 
            PrmParserStructures::isImproperSection(line) || PrmParserStructures::isNonbondedSection(line) || 
            line == "BOMLEV" || line == "WRNLEV" || line == "return") {
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        // Skip special directive lines like "HBOND CUTHB 0.5"
        if (tokens.size() >= 1 && (tokens[0] == "HBOND" || tokens[0] == "NBFIX")) {
            continue;
        }
        
        if (tokens.size() >= 4) {
            try {
                std::string type1 = tokens[0];
                std::string type2 = tokens[1];
                double epsilon = PrmParserStructures::safe_stod(tokens[2], "NBFIX epsilon for " + type1 + "-" + type2);
                double rmin = PrmParserStructures::safe_stod(tokens[3], "NBFIX Rmin for " + type1 + "-" + type2);
                
                // In CHARMM, NBFIX parameters are specified with full Rmin value
                // No need to multiply by 2 since we store the full Rmin value
                ff.add_nbfix(type1, type2, epsilon, rmin);
                
                if (PRMParser::debug_output) std::cerr << "Stored NBFIX for " << type1 << "-" << type2 
                    << ": epsilon = " << epsilon << ", Rmin = " << rmin << std::endl;
            } catch (const std::exception& e) {
                if (PRMParser::debug_output) std::cerr << "Warning: Skipping NBFIX line due to parsing error: " << line << std::endl;
            }
        }
    }
}

} // namespace io
} // namespace pygcmc