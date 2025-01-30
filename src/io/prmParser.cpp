// src/io/prmParser.cpp

#include "prmParser.hpp"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>

namespace pygcmc {

namespace {
    double safe_stod(const std::string& str, const std::string& context) {
        try {
            return std::stod(str);
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to convert '" + str + "' to double in " + context);
        }
    }

    int safe_stoi(const std::string& str, const std::string& context) {
        try {
            return std::stoi(str);
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to convert '" + str + "' to integer in " + context);
        }
    }
}

void PrmParser::parse_string(const std::string& content, ForceField& ff) {
    std::istringstream iss(content);
    PrmParser parser;
    parser.parseStream(iss, ff);
}

void PrmParser::parse_file(const std::string& filename, ForceField& ff) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open parameter file: " + filename);
    }
    PrmParser parser;
    parser.parseStream(file, ff);
}

void PrmParser::parse(const std::string& filename, ForceField& ff) {
    parse_file(filename, ff);
}

void PrmParser::parseStream(std::istream& input, ForceField& ff) {
    std::string line;
    bool inSection = false;
    std::string currentSection;
    
    while (std::getline(input, line)) {
        std::cerr << "Raw line: [" << line << "]" << std::endl;
        
        if (isCommentLine(line)) {
            std::cerr << "Skipping comment line" << std::endl;
            continue;
        }
        
        std::string cleanLine = removeComments(line);
        cleanLine = trim(cleanLine);
        std::cerr << "Cleaned line: [" << cleanLine << "]" << std::endl;
        
        if (cleanLine.empty()) continue;
        
        if (cleanLine == "END") {
            inSection = false;
            currentSection.clear();
            continue;
        }
        
        if (isAtomsSection(cleanLine)) {
            std::cerr << "Found ATOMS section" << std::endl;
            inSection = true;
            currentSection = "ATOMS";
            parseAtomsSection(input, ff);
        } else if (isBondsSection(cleanLine)) {
            std::cerr << "Found BONDS section" << std::endl;
            inSection = true;
            currentSection = "BONDS";
            parseBondsSection(input, ff);
        } else if (isAnglesSection(cleanLine)) {
            std::cerr << "Found ANGLES section" << std::endl;
            inSection = true;
            currentSection = "ANGLES";
            parseAnglesSection(input, ff);
        } else if (isDihedralsSection(cleanLine)) {
            std::cerr << "Found DIHEDRALS section" << std::endl;
            inSection = true;
            currentSection = "DIHEDRALS";
            parseDihedralsSection(input, ff);
        } else if (isImproperSection(cleanLine)) {
            std::cerr << "Found IMPROPER section" << std::endl;
            inSection = true;
            currentSection = "IMPROPER";
            parseImproperSection(input, ff);
        } else if (isNonbondedSection(cleanLine)) {
            std::cerr << "Found NONBONDED section" << std::endl;
            inSection = true;
            currentSection = "NONBONDED";
            parseNonbondedSection(input, ff, cleanLine);
        } else if (isNBFixSection(cleanLine)) {
            std::cerr << "Found NBFIX section" << std::endl;
            inSection = true;
            currentSection = "NBFIX";
            parseNBFixSection(input, ff);
        } else if (!inSection) {
            // Handle non-section content if needed
            auto tokens = tokenize(cleanLine);
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

void PrmParser::skipComments(std::istream& input) {
    while (input.peek() == '!' || input.peek() == '*' || input.peek() == '#' || input.peek() == '\n') {
        std::string line;
        std::getline(input, line);
    }
}

std::vector<std::string> PrmParser::tokenize(const std::string& line) {
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

bool PrmParser::isCommentLine(const std::string& line) {
    return line.empty() || line[0] == '!' || line[0] == '*' || line[0] == '#';
}

std::string PrmParser::removeComments(const std::string& line) {
    size_t commentPos = line.find_first_of("!*#");
    if (commentPos != std::string::npos) {
        return line.substr(0, commentPos);
    }
    return line;
}

std::string PrmParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::string PrmParser::readContinuationLine(std::istream& input, std::string firstLine) {
    std::cerr << "Reading continuation line starting with: [" << firstLine << "]" << std::endl;
    std::string fullLine;
    std::string currentLine = firstLine;
    
    while (true) {
        // Remove comments and trim the current line
        currentLine = removeComments(currentLine);
        currentLine = trim(currentLine);
        
        if (currentLine.empty()) {
            if (fullLine.empty()) {
                fullLine = currentLine;
            }
            break;
        }
        
        // Check if line ends with continuation character
        bool hasContinuation = false;
        if (!currentLine.empty() && currentLine.back() == '-') {
            hasContinuation = true;
            currentLine.pop_back();  // Remove the continuation character
            currentLine = trim(currentLine);  // Trim again after removing '-'
        }
        
        // Add the current line to the full line
        if (!fullLine.empty() && !currentLine.empty()) {
            fullLine += " ";  // Add space between continued lines
        }
        fullLine += currentLine;
        
        std::cerr << "Current full line: [" << fullLine << "]" << std::endl;
        
        // If no continuation character, we're done
        if (!hasContinuation) {
            break;
        }
        
        // Read next line
        if (!std::getline(input, currentLine)) {
            break;  // End of file
        }
        std::cerr << "Read continuation line: [" << currentLine << "]" << std::endl;
        
        // Skip comment lines in continuation
        while (isCommentLine(currentLine)) {
            if (!std::getline(input, currentLine)) {
                return fullLine;  // End of file
            }
            std::cerr << "Skipping comment in continuation: [" << currentLine << "]" << std::endl;
        }
    }
    
    std::cerr << "Final combined line: [" << fullLine << "]" << std::endl;
    return fullLine;
}

bool PrmParser::isAtomsSection(const std::string& line) {
    return line.find("ATOMS") != std::string::npos || line.find("MASS") != std::string::npos;
}

bool PrmParser::isBondsSection(const std::string& line) {
    return line.find("BONDS") != std::string::npos;
}

bool PrmParser::isAnglesSection(const std::string& line) {
    return line.find("ANGLES") != std::string::npos;
}

bool PrmParser::isDihedralsSection(const std::string& line) {
    return line.find("DIHEDRALS") != std::string::npos;
}

bool PrmParser::isImproperSection(const std::string& line) {
    return line.find("IMPROPER") != std::string::npos;
}

bool PrmParser::isNonbondedSection(const std::string& line) {
    return line.find("NONBONDED") != std::string::npos || 
           line.find("cutnb") != std::string::npos;
}

bool PrmParser::isNBFixSection(const std::string& line) {
    return line.find("NBFIX") != std::string::npos;
}

void PrmParser::parseAtomsSection(std::istream& input, ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || isBondsSection(line) || isAnglesSection(line) || 
            isDihedralsSection(line) || isImproperSection(line) || 
            isNonbondedSection(line) || isNBFixSection(line)) {
            break;
        }
        
        auto tokens = tokenize(line);
        // Only process lines that start with MASS and have the correct number of tokens
        if (tokens.size() >= 4 && tokens[0] == "MASS") {
            // Format: MASS -1 type mass [comment]
            std::string atomType = tokens[2];
            double mass = safe_stod(tokens[3], "atom mass for type " + atomType);
            ff.atom_masses[atomType] = mass;
        }
    }
}

void PrmParser::parseBondsSection(std::istream& input, ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || isAtomsSection(line) || isAnglesSection(line) || 
            isDihedralsSection(line) || isImproperSection(line) || 
            isNonbondedSection(line) || isNBFixSection(line)) {
            break;
        }
        
        auto tokens = tokenize(line);
        if (tokens.size() >= 4) {
            // Format: type1 type2 Kb b0
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            double kb = safe_stod(tokens[2], "bond Kb for " + type1 + "-" + type2);
            double b0 = safe_stod(tokens[3], "bond b0 for " + type1 + "-" + type2);
            
            auto key = make_type_pair(type1, type2);
            BondParams params{kb, b0};
            ff.bond_params[key] = params;
            // Add symmetric pair
            key = make_type_pair(type2, type1);
            ff.bond_params[key] = params;
        }
    }
}

void PrmParser::parseAnglesSection(std::istream& input, ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || isAtomsSection(line) || isBondsSection(line) || 
            isDihedralsSection(line) || isImproperSection(line) || 
            isNonbondedSection(line) || isNBFixSection(line)) {
            break;
        }
        
        auto tokens = tokenize(line);
        if (tokens.size() >= 5) {
            // Format: type1 type2 type3 Ktheta Theta0 [Kub S0]
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            std::string type3 = tokens[2];
            double ktheta = safe_stod(tokens[3], "angle Ktheta for " + type1 + "-" + type2 + "-" + type3);
            double theta0 = safe_stod(tokens[4], "angle theta0 for " + type1 + "-" + type2 + "-" + type3);
            
            double kub = 0.0;
            double s0 = 0.0;
            if (tokens.size() >= 7) {
                kub = safe_stod(tokens[5], "angle Kub for " + type1 + "-" + type2 + "-" + type3);
                s0 = safe_stod(tokens[6], "angle S0 for " + type1 + "-" + type2 + "-" + type3);
            }
            
            auto key = make_type_triple(type1, type2, type3);
            AngleParams params{ktheta, theta0, kub, s0};
            ff.angle_params[key] = params;
            // Add symmetric triple
            key = make_type_triple(type3, type2, type1);
            ff.angle_params[key] = params;
        }
    }
}

void PrmParser::parseDihedralsSection(std::istream& input, ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || isAtomsSection(line) || isBondsSection(line) || 
            isAnglesSection(line) || isImproperSection(line) || 
            isNonbondedSection(line) || isNBFixSection(line)) {
            break;
        }
        
        auto tokens = tokenize(line);
        if (tokens.size() >= 7) {
            // Format: type1 type2 type3 type4 Kchi n delta
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            std::string type3 = tokens[2];
            std::string type4 = tokens[3];
            double kchi = safe_stod(tokens[4], "dihedral Kchi for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            int n = safe_stoi(tokens[5], "dihedral n for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            double delta = safe_stod(tokens[6], "dihedral delta for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            
            auto key = make_type_quad(type1, type2, type3, type4);
            DihedralParams params{kchi, n, delta};
            ff.dihedral_params[key].push_back(params);
            // Add symmetric quad
            key = make_type_quad(type4, type3, type2, type1);
            ff.dihedral_params[key].push_back(params);
        }
    }
}

void PrmParser::parseImproperSection(std::istream& input, ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || isAtomsSection(line) || isBondsSection(line) || 
            isAnglesSection(line) || isDihedralsSection(line) || 
            isNonbondedSection(line) || isNBFixSection(line)) {
            break;
        }
        
        auto tokens = tokenize(line);
        if (tokens.size() >= 6) {
            // Format: type1 type2 type3 type4 Kpsi psi0
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            std::string type3 = tokens[2];
            std::string type4 = tokens[3];
            double kpsi = safe_stod(tokens[4], "improper Kpsi for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            double psi0 = safe_stod(tokens[5], "improper psi0 for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            
            auto key = make_type_quad(type1, type2, type3, type4);
            ImproperParams params{kpsi, psi0};
            ff.improper_params[key] = params;
        }
    }
}

void PrmParser::parseNBFixSection(std::istream& input, ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section end or new section
        if (line == "END" || isAtomsSection(line) || isBondsSection(line) || 
            isAnglesSection(line) || isDihedralsSection(line) || 
            isImproperSection(line) || isNonbondedSection(line) || 
            line == "BOMLEV" || line == "WRNLEV" || line == "return") {
            break;
        }
        
        auto tokens = tokenize(line);
        if (tokens.size() >= 4) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            try {
                double epsilon = safe_stod(tokens[2], "NBFIX epsilon for " + type1 + "-" + type2);
                // Store both directions to ensure symmetric lookup
                auto key = make_type_pair(type1, type2);
                ff.nbfix[key] = epsilon;
                
                std::cerr << "Stored NBFIX for " << type1 << "-" << type2 << ": epsilon = " << epsilon << std::endl;
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to parse NBFIX line: " + line + "\nError: " + e.what());
            }
        } else if (!tokens.empty()) {  // Only throw if line has some tokens but not enough
            throw std::runtime_error("Invalid NBFIX format in line: " + line);
        }
    }
}

void PrmParser::parseNonbondedSection(std::istream& input, ForceField& ff, const std::string& firstLine) {
    std::string line = firstLine;
    std::string fullLine = readContinuationLine(input, line);
    
    std::cerr << "=== Entering NONBONDED section parsing ===" << std::endl;
    std::cerr << "First line: [" << firstLine << "]" << std::endl;
    
    // Parse header parameters
    auto tokens = tokenize(fullLine);
    if (!tokens.empty()) {
        std::cerr << "Initial tokens:";
        for (const auto& token : tokens) {
            std::cerr << " [" << token << "]";
        }
        std::cerr << std::endl;
        
        // Skip the NONBONDED keyword
        if (tokens[0] == "NONBONDED") {
            std::cerr << "Skipping NONBONDED keyword" << std::endl;
            tokens.erase(tokens.begin());
        }
        
        // Process parameters
        for (size_t i = 0; i < tokens.size(); ++i) {
            std::cerr << "Processing parameter token: [" << tokens[i] << "]" << std::endl;
            
            if (tokens[i] == "nbxmod" && i + 1 < tokens.size()) {
                ff.nonbonded_params.nbxmod = safe_stoi(tokens[++i], "nbxmod");
                std::cerr << "Set nbxmod = " << ff.nonbonded_params.nbxmod << std::endl;
            } else if (tokens[i] == "cutnb" && i + 1 < tokens.size()) {
                ff.nonbonded_params.cutnb = safe_stod(tokens[++i], "cutnb");
                std::cerr << "Set cutnb = " << ff.nonbonded_params.cutnb << std::endl;
            } else if (tokens[i] == "ctofnb" && i + 1 < tokens.size()) {
                ff.nonbonded_params.ctofnb = safe_stod(tokens[++i], "ctofnb");
                std::cerr << "Set ctofnb = " << ff.nonbonded_params.ctofnb << std::endl;
            } else if (tokens[i] == "ctonnb" && i + 1 < tokens.size()) {
                ff.nonbonded_params.ctonnb = safe_stod(tokens[++i], "ctonnb");
                std::cerr << "Set ctonnb = " << ff.nonbonded_params.ctonnb << std::endl;
            } else if (tokens[i] == "eps" && i + 1 < tokens.size()) {
                ff.nonbonded_params.eps = safe_stod(tokens[++i], "eps");
                std::cerr << "Set eps = " << ff.nonbonded_params.eps << std::endl;
            } else if (tokens[i] == "e14fac" && i + 1 < tokens.size()) {
                ff.nonbonded_params.e14fac = safe_stod(tokens[++i], "e14fac");
                std::cerr << "Set e14fac = " << ff.nonbonded_params.e14fac << std::endl;
            } else if (tokens[i] == "wmin" && i + 1 < tokens.size()) {
                ff.nonbonded_params.wmin = safe_stod(tokens[++i], "wmin");
                std::cerr << "Set wmin = " << ff.nonbonded_params.wmin << std::endl;
            } else if (tokens[i] == "cdiel") {
                ff.nonbonded_params.cdiel = true;
                std::cerr << "Set cdiel = true" << std::endl;
            } else if (tokens[i] == "fshift") {
                ff.nonbonded_params.fshift = true;
                std::cerr << "Set fshift = true" << std::endl;
            } else if (tokens[i] == "vatom") {
                ff.nonbonded_params.vatom = true;
                std::cerr << "Set vatom = true" << std::endl;
            } else if (tokens[i] == "vdistance") {
                ff.nonbonded_params.vdistance = true;
                std::cerr << "Set vdistance = true" << std::endl;
            } else if (tokens[i] == "vfswitch") {
                ff.nonbonded_params.vfswitch = true;
                std::cerr << "Set vfswitch = true" << std::endl;
            }
        }
    }
    
    std::cerr << "\n=== Starting atom type parameters parsing ===" << std::endl;
    std::cerr << "Current lj_params map size: " << ff.lj_params.size() << std::endl;
    
    // Parse atom type parameters
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) {
            std::cerr << "Skipping empty line" << std::endl;
            continue;
        }
        
        std::cerr << "Processing cleaned line: [" << line << "]" << std::endl;
        tokens = tokenize(line);
        std::cerr << "Tokens:";
        for (const auto& token : tokens) {
            std::cerr << " [" << token << "]";
        }
        std::cerr << std::endl;
        
        if (tokens.empty()) continue;
        
        // Check for section end or new section
        if (line == "END" || isAtomsSection(line) || isBondsSection(line) || 
            isAnglesSection(line) || isDihedralsSection(line) || 
            isImproperSection(line)) {
            std::cerr << "Found section end marker: " << tokens[0] << std::endl;
            break;
        }
        
        // If we find NBFIX, parse it as a new section
        if (isNBFixSection(line)) {
            std::cerr << "Found NBFIX section" << std::endl;
            parseNBFixSection(input, ff);
            break;
        }
        
        // Format: atomType ignored epsilon Rmin [ignored ignored ignored ignored]
        if (tokens.size() >= 4) {
            std::string atomType = tokens[0];
            // Skip the ignored value (usually 0.0)
            try {
                double epsilon = safe_stod(tokens[2], "LJ epsilon for " + atomType);
                double rmin = safe_stod(tokens[3], "LJ Rmin for " + atomType);
                
                std::cerr << "\n*** Parsing atom type: " << atomType << " ***" << std::endl;
                std::cerr << "  epsilon = " << epsilon << std::endl;
                std::cerr << "  rmin = " << rmin << std::endl;
                
                // Store the parameters
                LJParams params{epsilon, rmin};
                ff.lj_params[atomType] = params;
                
                std::cerr << "Successfully stored " << atomType << " parameters:" << std::endl;
                std::cerr << "  Stored epsilon = " << ff.lj_params[atomType].epsilon << std::endl;
                std::cerr << "  Stored rmin = " << ff.lj_params[atomType].rmin << std::endl;
                std::cerr << "Current lj_params map size: " << ff.lj_params.size() << std::endl;
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to parse LJ parameters for " + atomType + ": " + e.what());
            }
        } else if (!tokens.empty() && !isCommentLine(line)) {
            // If we have tokens but not enough, and it's not a comment line, throw ValueError
            throw std::runtime_error("Malformed NONBONDED parameters in line: " + line + 
                                   "\nExpected at least 4 tokens, got " + std::to_string(tokens.size()));
        }
    }
    
    std::cerr << "\n=== Finished NONBONDED section parsing ===" << std::endl;
    std::cerr << "Final lj_params map size: " << ff.lj_params.size() << std::endl;
}

// Helper functions for parameter processing
std::pair<std::string, std::string> PrmParser::make_type_pair(
    const std::string& type1, const std::string& type2) const {
    return type1 < type2 ? 
        std::make_pair(type1, type2) : 
        std::make_pair(type2, type1);
}

std::tuple<std::string, std::string, std::string> PrmParser::make_type_triple(
    const std::string& type1, const std::string& type2, const std::string& type3) const {
    return std::make_tuple(type1, type2, type3);
}

std::tuple<std::string, std::string, std::string, std::string> PrmParser::make_type_quad(
    const std::string& type1, const std::string& type2, 
    const std::string& type3, const std::string& type4) const {
    return std::make_tuple(type1, type2, type3, type4);
}

} // namespace pygcmc

