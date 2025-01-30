// src/io/prmParser.cpp

#include "prmParser.hpp"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>

namespace pygcmc {

// Initialize static debug flag
bool PRMParser::debug_output = false;

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

void PRMParser::parse_string(const std::string& content, ForceField& ff) {
    std::istringstream iss(content);
    PRMParser parser;
    parser.parseStream(iss, ff);
}

void PRMParser::parse_file_to_forcefield(const std::string& filename, ForceField& ff) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open parameter file: " + filename);
    }
    PRMParser parser;
    parser.parseStream(file, ff);
}

ForceField PRMParser::parse_file(const std::string& filename) {
    ForceField ff;
    parse_file_to_forcefield(filename, ff);
    return ff;
}

void PRMParser::parse(const std::string& filename, ForceField& ff) {
    parse_file_to_forcefield(filename, ff);
}

void PRMParser::parseStream(std::istream& input, ForceField& ff) {
    std::string line;
    bool inSection = false;
    std::string currentSection;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw line: [" << line << "]" << std::endl;
        
        if (isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line" << std::endl;
            continue;
        }
        
        std::string cleanLine = removeComments(line);
        cleanLine = trim(cleanLine);
        if (debug_output) std::cerr << "Cleaned line: [" << cleanLine << "]" << std::endl;
        
        if (cleanLine.empty()) continue;
        
        if (cleanLine == "END") {
            inSection = false;
            currentSection.clear();
            continue;
        }
        
        if (isAtomsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ATOMS section" << std::endl;
            inSection = true;
            currentSection = "ATOMS";
            parseAtomsSection(input, ff);
        } else if (isBondsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found BONDS section" << std::endl;
            inSection = true;
            currentSection = "BONDS";
            parseBondsSection(input, ff);
        } else if (isAnglesSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ANGLES section" << std::endl;
            inSection = true;
            currentSection = "ANGLES";
            parseAnglesSection(input, ff);
        } else if (isDihedralsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found DIHEDRALS section" << std::endl;
            inSection = true;
            currentSection = "DIHEDRALS";
            parseDihedralsSection(input, ff);
        } else if (isImproperSection(cleanLine)) {
            if (debug_output) std::cerr << "Found IMPROPER section" << std::endl;
            inSection = true;
            currentSection = "IMPROPER";
            parseImproperSection(input, ff);
        } else if (isNonbondedSection(cleanLine)) {
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
                    if (isCommentLine(prevLine)) continue;
                    prevLine = removeComments(prevLine);
                    prevLine = trim(prevLine);
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
        } else if (isNBFixSection(cleanLine)) {
            if (debug_output) std::cerr << "Found NBFIX section" << std::endl;
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

void PRMParser::skipComments(std::istream& input) {
    while (input.peek() == '!' || input.peek() == '*' || input.peek() == '#' || input.peek() == '\n') {
        std::string line;
        std::getline(input, line);
    }
}

std::vector<std::string> PRMParser::tokenize(const std::string& line) {
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

bool PRMParser::isCommentLine(const std::string& line) {
    return line.empty() || line[0] == '!' || line[0] == '*' || line[0] == '#';
}

std::string PRMParser::removeComments(const std::string& line) {
    size_t commentPos = line.find_first_of("!*#");
    if (commentPos != std::string::npos) {
        return line.substr(0, commentPos);
    }
    return line;
}

std::string PRMParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::string PRMParser::readContinuationLine(std::istream& input, std::string firstLine) {
    if (debug_output) std::cerr << "Reading continuation line starting with: [" << firstLine << "]" << std::endl;
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
        if (debug_output) std::cerr << "Read continuation line: [" << currentLine << "]" << std::endl;
        
        while (isCommentLine(currentLine)) {
            if (!std::getline(input, currentLine)) {
                return fullLine;
            }
            if (debug_output) std::cerr << "Skipping comment in continuation: [" << currentLine << "]" << std::endl;
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
        
        if (debug_output) std::cerr << "Current full line: [" << fullLine << "]" << std::endl;
    }
    
    if (debug_output) std::cerr << "Final combined line: [" << fullLine << "]" << std::endl;
    return fullLine;
}

bool PRMParser::isAtomsSection(const std::string& line) {
    return line.find("ATOMS") != std::string::npos || line.find("MASS") != std::string::npos;
}

bool PRMParser::isBondsSection(const std::string& line) {
    bool result = line.find("BONDS") != std::string::npos;
    if (result && debug_output) {
        std::cerr << "Found BONDS section marker: [" << line << "]" << std::endl;
    }
    return result;
}

bool PRMParser::isAnglesSection(const std::string& line) {
    return line.find("ANGLES") != std::string::npos;
}

bool PRMParser::isDihedralsSection(const std::string& line) {
    return line.find("DIHEDRALS") != std::string::npos;
}

bool PRMParser::isImproperSection(const std::string& line) {
    return line.find("IMPROPER") != std::string::npos;
}

bool PRMParser::isNonbondedSection(const std::string& line) {
    return line.find("NONBONDED") != std::string::npos || 
           line.find("cutnb") != std::string::npos;
}

bool PRMParser::isNBFixSection(const std::string& line) {
    return line.find("NBFIX") != std::string::npos;
}

void PRMParser::parseAtomsSection(std::istream& input, ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering ATOMS/MASS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw atom line: [" << line << "]" << std::endl;
        
        if (isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in atoms section" << std::endl;
            continue;
        }
        
        std::string fullLine = readContinuationLine(input, line);
        fullLine = removeComments(fullLine);
        fullLine = trim(fullLine);
        
        if (fullLine.empty()) {
            if (debug_output) std::cerr << "Skipping empty line in atoms section" << std::endl;
            continue;
        }
        
        if (debug_output) std::cerr << "Processed line: [" << fullLine << "]" << std::endl;
        
        // Check for section end
        if (fullLine == "END" || isBondsSection(fullLine) || isAnglesSection(fullLine) || 
            isDihedralsSection(fullLine) || isImproperSection(fullLine) || 
            isNonbondedSection(fullLine) || isNBFixSection(fullLine)) {
            if (debug_output) std::cerr << "Found section end marker: " << fullLine << std::endl;
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = tokenize(fullLine);
        if (debug_output) {
            std::cerr << "Tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.size() >= 4 && tokens[0] == "MASS") {
            std::string atomType = tokens[2];
            double mass = safe_stod(tokens[3], "atom mass for type " + atomType);
            ff.atom_masses[atomType] = mass;
            if (debug_output) std::cerr << "Added atom mass: " << atomType << " = " << mass << std::endl;
        } else {
            if (debug_output) std::cerr << "Skipping line: not a valid MASS entry" << std::endl;
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished ATOMS/MASS section parsing ===" << std::endl;
        std::cerr << "Final atom_masses size: " << ff.atom_masses.size() << std::endl;
        std::cerr << "\nStored atom masses:" << std::endl;
        for (const auto& pair : ff.atom_masses) {
            std::cerr << pair.first << " = " << pair.second << std::endl;
        }
    }
}

void PRMParser::parseBondsSection(std::istream& input, ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering BONDS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw bond line: [" << line << "]" << std::endl;
        
        if (isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in bonds section" << std::endl;
            continue;
        }
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) {
            if (debug_output) std::cerr << "Skipping empty line in bonds section" << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || isAtomsSection(line) || isAnglesSection(line) || 
            isDihedralsSection(line) || isImproperSection(line) || 
            isNonbondedSection(line) || isNBFixSection(line)) {
            if (debug_output) std::cerr << "Found section end marker in bonds: " << line << std::endl;
            // Put the line back so it can be read by the next section parser
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = tokenize(line);
        if (debug_output) {
            std::cerr << "Bond tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.size() >= 4) {
            // Format: type1 type2 Kb b0
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            double kb = safe_stod(tokens[2], "bond Kb for " + type1 + "-" + type2);
            double b0 = safe_stod(tokens[3], "bond b0 for " + type1 + "-" + type2);
            
            if (debug_output) {
                std::cerr << "\nProcessing bond: " << type1 << "-" << type2 << std::endl;
                std::cerr << "  kb = " << kb << ", b0 = " << b0 << std::endl;
            }
            
            // Store bond parameters in a consistent order
            auto key = ForceField::makeTypePair(type1, type2);
            if (debug_output) {
                std::cerr << "  Storing with key: (" << key.first << ", " << key.second << ")" << std::endl;
            }
            
            BondParams params{kb, b0};
            ff.bond_params[key] = params;
            
            if (debug_output) {
                std::cerr << "  Current bond_params size: " << ff.bond_params.size() << std::endl;
            }
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished BONDS section parsing ===" << std::endl;
        std::cerr << "Final bond_params size: " << ff.bond_params.size() << std::endl;
    }
}

void PRMParser::parseAnglesSection(std::istream& input, ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering ANGLES section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw angle line: [" << line << "]" << std::endl;
        
        if (isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in angles section" << std::endl;
            continue;
        }
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        // Check for section end
        if (line == "END" || isAtomsSection(line) || isBondsSection(line) || 
            isDihedralsSection(line) || isImproperSection(line) || 
            isNonbondedSection(line) || isNBFixSection(line)) {
            // Put the line back so it can be read by the next section parser
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = tokenize(line);
        if (tokens.size() >= 5) {
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
            
            if (debug_output) {
                std::cerr << "\nProcessing angle: " << type1 << "-" << type2 << "-" << type3 << std::endl;
                std::cerr << "  ktheta = " << ktheta << ", theta0 = " << theta0 << std::endl;
            }
            
            AngleParams params{ktheta, theta0, kub, s0};
            
            // Store all possible permutations
            std::vector<std::tuple<std::string, std::string, std::string>> keys = {
                make_type_triple(type1, type2, type3),
                make_type_triple(type3, type2, type1),
                std::make_tuple(type2, type1, type3),
                std::make_tuple(type2, type3, type1)
            };
            
            for (const auto& key : keys) {
                ff.angle_params[key] = params;
                if (debug_output) {
                    std::cerr << "  Storing with key: (" 
                             << std::get<0>(key) << ", "
                             << std::get<1>(key) << ", "
                             << std::get<2>(key) << ")" << std::endl;
                }
            }
            
            if (debug_output) {
                std::cerr << "  Current angle_params size: " << ff.angle_params.size() << std::endl;
            }
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished ANGLES section parsing ===" << std::endl;
        std::cerr << "Final angle_params size: " << ff.angle_params.size() << std::endl;
        std::cerr << "\nStored angle parameter keys:" << std::endl;
        for (const auto& pair : ff.angle_params) {
            const auto& key = pair.first;
            std::cerr << "(" << std::get<0>(key) << ", "
                     << std::get<1>(key) << ", "
                     << std::get<2>(key) << ")" << std::endl;
        }
    }
}

void PRMParser::parseDihedralsSection(std::istream& input, ForceField& ff) {
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
            // Put the line back so it can be read by the next section parser
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
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

void PRMParser::parseImproperSection(std::istream& input, ForceField& ff) {
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
            // Put the line back so it can be read by the next section parser
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = tokenize(line);
        if (tokens.size() >= 6) {
            try {
                // Format: type1 type2 type3 type4 Kpsi psi0
                std::string type1 = tokens[0];
                std::string type2 = tokens[1];
                std::string type3 = tokens[2];
                std::string type4 = tokens[3];
                // Only try to convert the last two tokens to numbers
                double kpsi = safe_stod(tokens[4], "improper Kpsi");
                double psi0 = safe_stod(tokens[5], "improper psi0");
                
                auto key = make_type_quad(type1, type2, type3, type4);
                ImproperParams params{kpsi, psi0};
                ff.improper_params[key] = params;
            } catch (const std::exception& e) {
                // Skip lines that can't be parsed properly
                if (debug_output) std::cerr << "Warning: Skipping improper line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

void PRMParser::parseNBFixSection(std::istream& input, ForceField& ff) {
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
            // Put the line back so it can be read by the next section parser
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = tokenize(line);
        // Skip special directive lines like "HBOND CUTHB 0.5"
        if (tokens.size() >= 1 && (tokens[0] == "HBOND" || tokens[0] == "NBFIX")) {
            continue;
        }
        
        if (tokens.size() >= 4) {
            try {
                std::string type1 = tokens[0];
                std::string type2 = tokens[1];
                double epsilon = safe_stod(tokens[2], "NBFIX epsilon for " + type1 + "-" + type2);
                // Store both directions to ensure symmetric lookup
                auto key = make_type_pair(type1, type2);
                ff.nbfix[key] = epsilon;
                key = make_type_pair(type2, type1);
                ff.nbfix[key] = epsilon;
                
                if (debug_output) std::cerr << "Stored NBFIX for " << type1 << "-" << type2 << ": epsilon = " << epsilon << std::endl;
            } catch (const std::exception& e) {
                // Skip lines that can't be parsed properly
                if (debug_output) std::cerr << "Warning: Skipping NBFIX line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

void PRMParser::parseNonbondedSection(std::istream& input, ForceField& ff, const std::string& firstLine) {
    std::string line = firstLine;
    std::string fullLine = readContinuationLine(input, line);
    
    if (debug_output) std::cerr << "=== Entering NONBONDED section parsing ===" << std::endl;
    if (debug_output) std::cerr << "First line: [" << firstLine << "]" << std::endl;
    
    // Parse header parameters
    auto tokens = tokenize(fullLine);
    if (!tokens.empty()) {
        if (debug_output) {
            std::cerr << "Initial tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        // Skip the NONBONDED keyword
        if (tokens[0] == "NONBONDED") {
            if (debug_output) std::cerr << "Skipping NONBONDED keyword" << std::endl;
            tokens.erase(tokens.begin());
        }
        
        // Process parameters
        for (size_t i = 0; i < tokens.size(); ++i) {
            if (debug_output) std::cerr << "Processing parameter token: [" << tokens[i] << "]" << std::endl;
            
            if (tokens[i] == "nbxmod" && i + 1 < tokens.size()) {
                ff.nonbonded_params.nbxmod = safe_stoi(tokens[++i], "nbxmod");
                if (debug_output) std::cerr << "Set nbxmod = " << ff.nonbonded_params.nbxmod << std::endl;
            } else if (tokens[i] == "cutnb" && i + 1 < tokens.size()) {
                ff.nonbonded_params.cutnb = safe_stod(tokens[++i], "cutnb");
                if (debug_output) std::cerr << "Set cutnb = " << ff.nonbonded_params.cutnb << std::endl;
            } else if (tokens[i] == "ctofnb" && i + 1 < tokens.size()) {
                ff.nonbonded_params.ctofnb = safe_stod(tokens[++i], "ctofnb");
                if (debug_output) std::cerr << "Set ctofnb = " << ff.nonbonded_params.ctofnb << std::endl;
            } else if (tokens[i] == "ctonnb" && i + 1 < tokens.size()) {
                ff.nonbonded_params.ctonnb = safe_stod(tokens[++i], "ctonnb");
                if (debug_output) std::cerr << "Set ctonnb = " << ff.nonbonded_params.ctonnb << std::endl;
            } else if (tokens[i] == "eps" && i + 1 < tokens.size()) {
                ff.nonbonded_params.eps = safe_stod(tokens[++i], "eps");
                if (debug_output) std::cerr << "Set eps = " << ff.nonbonded_params.eps << std::endl;
            } else if (tokens[i] == "e14fac" && i + 1 < tokens.size()) {
                ff.nonbonded_params.e14fac = safe_stod(tokens[++i], "e14fac");
                if (debug_output) std::cerr << "Set e14fac = " << ff.nonbonded_params.e14fac << std::endl;
            } else if (tokens[i] == "wmin" && i + 1 < tokens.size()) {
                ff.nonbonded_params.wmin = safe_stod(tokens[++i], "wmin");
                if (debug_output) std::cerr << "Set wmin = " << ff.nonbonded_params.wmin << std::endl;
            } else if (tokens[i] == "cdiel") {
                ff.nonbonded_params.cdiel = true;
                if (debug_output) std::cerr << "Set cdiel = true" << std::endl;
            } else if (tokens[i] == "fshift") {
                ff.nonbonded_params.fshift = true;
                if (debug_output) std::cerr << "Set fshift = true" << std::endl;
            } else if (tokens[i] == "vatom") {
                ff.nonbonded_params.vatom = true;
                if (debug_output) std::cerr << "Set vatom = true" << std::endl;
            } else if (tokens[i] == "vdistance") {
                ff.nonbonded_params.vdistance = true;
                if (debug_output) std::cerr << "Set vdistance = true" << std::endl;
            } else if (tokens[i] == "vfswitch") {
                ff.nonbonded_params.vfswitch = true;
                if (debug_output) std::cerr << "Set vfswitch = true" << std::endl;
            }
        }
    }
    
    if (debug_output) std::cerr << "\n=== Starting atom type parameters parsing ===" << std::endl;
    if (debug_output) std::cerr << "Current lj_params map size: " << ff.lj_params.size() << std::endl;
    
    // Parse atom type parameters
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) {
            if (debug_output) std::cerr << "Skipping empty line" << std::endl;
            continue;
        }
        
        if (debug_output) std::cerr << "Processing cleaned line: [" << line << "]" << std::endl;
        tokens = tokenize(line);
        if (debug_output) {
            std::cerr << "Tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.empty()) continue;
        
        // Check for section end or new section
        if (line == "END" || isAtomsSection(line) || isBondsSection(line) || 
            isAnglesSection(line) || isDihedralsSection(line) || 
            isImproperSection(line)) {
            if (debug_output) std::cerr << "Found section end marker: " << tokens[0] << std::endl;
            break;
        }
        
        // If we find NBFIX, parse it as a new section
        if (isNBFixSection(line)) {
            if (debug_output) std::cerr << "Found NBFIX section" << std::endl;
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
                
                if (debug_output) {
                    std::cerr << "\n*** Parsing atom type: " << atomType << " ***" << std::endl;
                    std::cerr << "  epsilon = " << epsilon << std::endl;
                    std::cerr << "  rmin = " << rmin << std::endl;
                }
                
                // Store the parameters
                LJParams params{epsilon, rmin};
                ff.lj_params[atomType] = params;
                
                if (debug_output) {
                    std::cerr << "Successfully stored " << atomType << " parameters:" << std::endl;
                    std::cerr << "  Stored epsilon = " << ff.lj_params[atomType].epsilon << std::endl;
                    std::cerr << "  Stored rmin = " << ff.lj_params[atomType].rmin << std::endl;
                    std::cerr << "Current lj_params map size: " << ff.lj_params.size() << std::endl;
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to parse LJ parameters for " + atomType + ": " + e.what());
            }
        } else if (!tokens.empty() && !isCommentLine(line)) {
            // If we have tokens but not enough, and it's not a comment line, throw ValueError
            throw std::runtime_error("Malformed NONBONDED parameters in line: " + line + 
                                   "\nExpected at least 4 tokens, got " + std::to_string(tokens.size()));
        }
    }
    
    if (debug_output) std::cerr << "\n=== Finished NONBONDED section parsing ===" << std::endl;
    if (debug_output) std::cerr << "Final lj_params map size: " << ff.lj_params.size() << std::endl;
}

// Helper functions for parameter processing
std::pair<std::string, std::string> PRMParser::make_type_pair(
    const std::string& type1, const std::string& type2) const {
    return type1 < type2 ? 
        std::make_pair(type1, type2) : 
        std::make_pair(type2, type1);
}

std::tuple<std::string, std::string, std::string> PRMParser::make_type_triple(
    const std::string& type1, const std::string& type2, const std::string& type3) const {
    // For angle parameters, we need to handle both symmetric and alternative representations
    // Create a vector of all possible representations
    std::vector<std::tuple<std::string, std::string, std::string>> keys = {
        std::make_tuple(type1, type2, type3),  // original order
        std::make_tuple(type3, type2, type1),  // symmetric order
        std::make_tuple(type2, type1, type3),  // alternative representation
        std::make_tuple(type2, type3, type1)   // symmetric alternative representation
    };
    
    // Return the lexicographically smallest key to ensure consistency
    return *std::min_element(keys.begin(), keys.end());
}

std::tuple<std::string, std::string, std::string, std::string> PRMParser::make_type_quad(
    const std::string& type1, const std::string& type2, 
    const std::string& type3, const std::string& type4) const {
    return std::make_tuple(type1, type2, type3, type4);
}

} // namespace pygcmc

