#include "prmParser.hpp"
#include <algorithm>
#include <cctype>
#include <iomanip>

namespace pygcmc {

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
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '!' || line[0] == '*' || line[0] == '#') {
            continue;
        }
        
        if (isNonbondedSection(line)) {
            // Include the NONBONDED line in parsing
            std::string fullLine = line;
            parseNonbondedSection(input, ff, fullLine);
        } else if (isNBFixSection(line)) {
            parseNBFixSection(input, ff);
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

bool PrmParser::isNonbondedSection(const std::string& line) {
    return line.find("NONBONDED") != std::string::npos;
}

bool PrmParser::isNBFixSection(const std::string& line) {
    return line.find("NBFIX") != std::string::npos;
}

void PrmParser::parseNonbondedSection(std::istream& input, ForceField& ff, const std::string& firstLine) {
    NonbondedParams& params = ff.get_nonbonded_params();
    // Initialize boolean parameters to false
    params.cdiel = false;
    params.fshift = false;
    params.vatom = false;
    params.vdistance = false;
    params.vfswitch = false;
    
    std::string fullLine = firstLine;  // Start with the NONBONDED line
    bool hasContinuation = false;
    
    // Check if first line has continuation
    if (!fullLine.empty()) {
        size_t commentPos = fullLine.find_first_of("!*#");
        if (commentPos != std::string::npos) {
            fullLine = fullLine.substr(0, commentPos);
        }
        
        // Trim right whitespace
        while (!fullLine.empty() && std::isspace(fullLine.back())) {
            fullLine.pop_back();
        }
        
        if (!fullLine.empty() && fullLine.back() == '-') {
            fullLine.pop_back();  // Remove continuation character
            hasContinuation = true;
        }
    }
    
    // Read continuation lines if any
    std::string line;
    while (hasContinuation && std::getline(input, line)) {
        if (line.empty() || line[0] == '!' || line[0] == '*' || line[0] == '#') {
            continue;
        }
        
        // Remove comments
        size_t commentPos = line.find_first_of("!*#");
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        
        // Trim right whitespace
        while (!line.empty() && std::isspace(line.back())) {
            line.pop_back();
        }
        
        hasContinuation = false;
        if (!line.empty() && line.back() == '-') {
            line.pop_back();  // Remove continuation character
            hasContinuation = true;
        }
        
        fullLine += " " + line;
    }
    
    auto tokens = tokenize(fullLine);
    
    // Skip the NONBONDED keyword if present
    size_t startIdx = 0;
    if (!tokens.empty() && tokens[0] == "NONBONDED") {
        startIdx = 1;
    }
    
    // Parse parameters
    for (size_t i = startIdx; i < tokens.size(); ++i) {
        const std::string& token = tokens[i];
        if (token == "nbxmod" && i + 1 < tokens.size()) 
            params.nbxmod = std::stoi(tokens[i+1]);
        else if (token == "cdiel") 
            params.cdiel = true;
        else if (token == "fshift") 
            params.fshift = true;
        else if (token == "vatom") 
            params.vatom = true;
        else if (token == "vdistance") 
            params.vdistance = true;
        else if (token == "vfswitch") 
            params.vfswitch = true;
        else if (token == "cutnb" && i + 1 < tokens.size()) 
            params.cutnb = std::stod(tokens[i+1]);
        else if (token == "ctofnb" && i + 1 < tokens.size()) 
            params.ctofnb = std::stod(tokens[i+1]);
        else if (token == "ctonnb" && i + 1 < tokens.size()) 
            params.ctonnb = std::stod(tokens[i+1]);
        else if (token == "eps" && i + 1 < tokens.size()) 
            params.eps = std::stod(tokens[i+1]);
        else if (token == "e14fac" && i + 1 < tokens.size()) 
            params.e14fac = std::stod(tokens[i+1]);
        else if (token == "wmin" && i + 1 < tokens.size()) 
            params.wmin = std::stod(tokens[i+1]);
    }
    
    // Parse atom type parameters
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '!' || line[0] == '*' || line[0] == '#') {
            continue;
        }
        
        tokens = tokenize(line);
        if (tokens.empty()) continue;
        
        if (tokens[0] == "END" || tokens[0] == "NBFIX") break;
        
        if (tokens.size() >= 4) {
            std::string atomType = tokens[0];
            LJParams ljParams;
            ljParams.epsilon = std::stod(tokens[2]);  // epsilon is in the third column
            ljParams.rmin = std::stod(tokens[3]);     // rmin is in the fourth column
            ff.add_lj_params(atomType, ljParams);
        }
    }
}

void PrmParser::parseNBFixSection(std::istream& input, ForceField& ff) {
    std::string line;
    
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '!' || line[0] == '*' || line[0] == '#') {
            continue;
        }
        
        auto tokens = tokenize(line);
        if (tokens.empty()) continue;
        
        if (tokens[0] == "END") break;
        
        if (tokens.size() >= 4) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            double epsilon = std::stod(tokens[2]);
            double rmin = std::stod(tokens[3]);
            ff.add_nbfix(type1, type2, epsilon, rmin);
        }
    }
}

} // namespace pygcmc
