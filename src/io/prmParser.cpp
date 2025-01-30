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
    bool inNonbondedSection = false;
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
        
        if (isAtomsSection(cleanLine)) {
            std::cerr << "Found ATOMS section" << std::endl;
            parseAtomsSection(input, ff);
        } else if (isBondsSection(cleanLine)) {
            std::cerr << "Found BONDS section" << std::endl;
            parseBondsSection(input, ff);
        } else if (isAnglesSection(cleanLine)) {
            std::cerr << "Found ANGLES section" << std::endl;
            parseAnglesSection(input, ff);
        } else if (isDihedralsSection(cleanLine)) {
            std::cerr << "Found DIHEDRALS section" << std::endl;
            parseDihedralsSection(input, ff);
        } else if (isImproperSection(cleanLine)) {
            std::cerr << "Found IMPROPER section" << std::endl;
            parseImproperSection(input, ff);
        } else if (isNonbondedSection(cleanLine)) {
            std::cerr << "Found NONBONDED section" << std::endl;
            parseNonbondedSection(input, ff, cleanLine);
        } else if (isNBFixSection(cleanLine)) {
            std::cerr << "Found NBFIX section" << std::endl;
            parseNBFixSection(input, ff);
        } else if (inNonbondedSection) {
            auto tokens = tokenize(cleanLine);
            if (!tokens.empty()) {
                if (tokens[0] == "NBFIX" || tokens[0] == "END") {
                    std::cerr << "End of NONBONDED section" << std::endl;
                    inNonbondedSection = false;
                } else if (tokens.size() >= 4) {
                    // 这里处理参数...
                }
            }
            continue;
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

void PrmParser::parseNonbondedSection(std::istream& input, ForceField& ff, const std::string& firstLine) {
    std::cerr << "\n=== Entering NONBONDED section parsing ===" << std::endl;
    std::cerr << "First line: [" << firstLine << "]" << std::endl;
    
    // Initialize boolean parameters to false
    ff.nonbonded_params.cdiel = false;
    ff.nonbonded_params.fshift = false;
    ff.nonbonded_params.vatom = false;
    ff.nonbonded_params.vdistance = false;
    ff.nonbonded_params.vfswitch = false;
    
    // Parse the first line which contains the NONBONDED parameters
    std::string fullLine = readContinuationLine(input, firstLine);
    std::cerr << "After continuation, full line: [" << fullLine << "]" << std::endl;
    
    auto tokens = tokenize(fullLine);
    std::cerr << "Initial tokens:";
    for (const auto& token : tokens) {
        std::cerr << " [" << token << "]";
    }
    std::cerr << std::endl;
    
    // Skip the NONBONDED keyword if present
    size_t startIdx = 0;
    if (!tokens.empty() && tokens[0] == "NONBONDED") {
        startIdx = 1;
        std::cerr << "Skipping NONBONDED keyword" << std::endl;
    }
    
    // Parse parameters
    for (size_t i = startIdx; i < tokens.size(); ++i) {
        const std::string& token = tokens[i];
        std::cerr << "Processing parameter token: [" << token << "]" << std::endl;
        if (token == "nbxmod" && i + 1 < tokens.size()) {
            ff.nonbonded_params.nbxmod = safe_stoi(tokens[i+1], "nbxmod");
            std::cerr << "Set nbxmod = " << ff.nonbonded_params.nbxmod << std::endl;
        }
        else if (token == "cdiel") {
            ff.nonbonded_params.cdiel = true;
            std::cerr << "Set cdiel = true" << std::endl;
        }
        else if (token == "fshift") {
            ff.nonbonded_params.fshift = true;
            std::cerr << "Set fshift = true" << std::endl;
        }
        else if (token == "vatom") {
            ff.nonbonded_params.vatom = true;
            std::cerr << "Set vatom = true" << std::endl;
        }
        else if (token == "vdistance") {
            ff.nonbonded_params.vdistance = true;
            std::cerr << "Set vdistance = true" << std::endl;
        }
        else if (token == "vfswitch") {
            ff.nonbonded_params.vfswitch = true;
            std::cerr << "Set vfswitch = true" << std::endl;
        }
        else if (token == "cutnb" && i + 1 < tokens.size()) {
            ff.nonbonded_params.cutnb = safe_stod(tokens[i+1], "cutnb");
            std::cerr << "Set cutnb = " << ff.nonbonded_params.cutnb << std::endl;
        }
        else if (token == "ctofnb" && i + 1 < tokens.size()) {
            ff.nonbonded_params.ctofnb = safe_stod(tokens[i+1], "ctofnb");
            std::cerr << "Set ctofnb = " << ff.nonbonded_params.ctofnb << std::endl;
        }
        else if (token == "ctonnb" && i + 1 < tokens.size()) {
            ff.nonbonded_params.ctonnb = safe_stod(tokens[i+1], "ctonnb");
            std::cerr << "Set ctonnb = " << ff.nonbonded_params.ctonnb << std::endl;
        }
        else if (token == "eps" && i + 1 < tokens.size()) {
            ff.nonbonded_params.eps = safe_stod(tokens[i+1], "eps");
            std::cerr << "Set eps = " << ff.nonbonded_params.eps << std::endl;
        }
        else if (token == "e14fac" && i + 1 < tokens.size()) {
            ff.nonbonded_params.e14fac = safe_stod(tokens[i+1], "e14fac");
            std::cerr << "Set e14fac = " << ff.nonbonded_params.e14fac << std::endl;
        }
        else if (token == "wmin" && i + 1 < tokens.size()) {
            ff.nonbonded_params.wmin = safe_stod(tokens[i+1], "wmin");
            std::cerr << "Set wmin = " << ff.nonbonded_params.wmin << std::endl;
        }
    }
    
    // Parse atom type parameters
    std::string line;
    std::cerr << "\n=== Starting atom type parameters parsing ===" << std::endl;
    std::cerr << "Current lj_params map size: " << ff.lj_params.size() << std::endl;
    
    while (std::getline(input, line)) {
        std::cerr << "\nReading line: [" << line << "]" << std::endl;
        if (isCommentLine(line)) {
            std::cerr << "Skipping comment line" << std::endl;
            continue;
        }
        
        // Handle continuation lines
        std::string originalLine = line;
        line = readContinuationLine(input, line);
        if (line != originalLine) {
            std::cerr << "After continuation processing: [" << line << "]" << std::endl;
        }
        
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
        
        // Check for section end
        if (tokens[0] == "END" || tokens[0] == "NBFIX" || 
            isAtomsSection(line) || isBondsSection(line) || 
            isAnglesSection(line) || isDihedralsSection(line) || 
            isImproperSection(line) || isNonbondedSection(line)) {
            std::cerr << "Found section end marker: " << tokens[0] << std::endl;
            break;
        }
        
        // Format: atomType ignored epsilon Rmin [ignored ignored ignored ignored]
        if (tokens.size() >= 4) {
            std::string atomType = tokens[0];
            // Skip the ignored value (usually 0.0)
            double epsilon = safe_stod(tokens[2], "LJ epsilon for " + atomType);
            double rmin = safe_stod(tokens[3], "LJ rmin for " + atomType);
            
            std::cerr << "\n*** Parsing atom type: " << atomType << " ***" << std::endl;
            std::cerr << "  epsilon = " << epsilon << std::endl;
            std::cerr << "  rmin = " << rmin << std::endl;
            
            LJParams ljParams;
            ljParams.epsilon = epsilon;
            ljParams.rmin = rmin;
            ff.lj_params[atomType] = ljParams;
            
            // Verify storage
            auto it = ff.lj_params.find(atomType);
            if (it != ff.lj_params.end()) {
                std::cerr << "Successfully stored " << atomType << " parameters:" << std::endl;
                std::cerr << "  Stored epsilon = " << it->second.epsilon << std::endl;
                std::cerr << "  Stored rmin = " << it->second.rmin << std::endl;
                std::cerr << "Current lj_params map size: " << ff.lj_params.size() << std::endl;
            } else {
                std::cerr << "WARNING: Failed to store " << atomType << " in map!" << std::endl;
            }
        } else {
            std::cerr << "Line has insufficient tokens (" << tokens.size() << " < 4), skipping" << std::endl;
        }
    }
    std::cerr << "\n=== Finished NONBONDED section parsing ===" << std::endl;
    std::cerr << "Final lj_params map size: " << ff.lj_params.size() << std::endl;
}

void PrmParser::parseNBFixSection(std::istream& input, ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (isCommentLine(line)) continue;
        
        line = removeComments(line);
        line = trim(line);
        if (line.empty()) continue;
        
        auto tokens = tokenize(line);
        if (tokens.empty()) continue;
        
        if (tokens[0] == "END") break;
        
        if (tokens.size() >= 4) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            double epsilon = safe_stod(tokens[2], "NBFIX epsilon for " + type1 + "-" + type2);
            // rmin is not used in current implementation
            
            auto key = make_type_pair(type1, type2);
            ff.nbfix[key] = epsilon;
            // Add symmetric pair
            key = make_type_pair(type2, type1);
            ff.nbfix[key] = epsilon;
        }
    }
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

