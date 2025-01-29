#include "prmParser.hpp"
#include <algorithm>
#include <cctype>

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
        skipComments(input);
        
        if (isNonbondedSection(line)) {
            parseNonbondedSection(input, ff);
        } else if (isNBFixSection(line)) {
            parseNBFixSection(input, ff);
        }
    }
}

void PrmParser::skipComments(std::istream& input) {
    std::string line;
    while (input.peek() == '!' || input.peek() == '*' || input.peek() == '#') {
        std::getline(input, line);
    }
}

std::vector<std::string> PrmParser::tokenize(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    
    while (iss >> token) {
        if (token[0] == '!' || token[0] == '#') break;  // Stop at comments
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

void PrmParser::parseNonbondedSection(std::istream& input, ForceField& ff) {
    NonbondedParams& params = ff.get_nonbonded_params();
    std::string line;
    
    // Parse header line
    std::getline(input, line);
    auto tokens = tokenize(line);
    
    // Parse parameters
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (tokens[i] == "nbxmod") params.nbxmod = std::stoi(tokens[i+1]);
        else if (tokens[i] == "cdiel") params.cdiel = true;
        else if (tokens[i] == "fshift") params.fshift = true;
        else if (tokens[i] == "vatom") params.vatom = true;
        else if (tokens[i] == "vdistance") params.vdistance = true;
        else if (tokens[i] == "vfswitch") params.vfswitch = true;
        else if (tokens[i] == "cutnb") params.cutnb = std::stod(tokens[i+1]);
        else if (tokens[i] == "ctofnb") params.ctofnb = std::stod(tokens[i+1]);
        else if (tokens[i] == "ctonnb") params.ctonnb = std::stod(tokens[i+1]);
        else if (tokens[i] == "eps") params.eps = std::stod(tokens[i+1]);
        else if (tokens[i] == "e14fac") params.e14fac = std::stod(tokens[i+1]);
        else if (tokens[i] == "wmin") params.wmin = std::stod(tokens[i+1]);
    }
    
    // Parse atom type parameters
    while (std::getline(input, line)) {
        skipComments(input);
        tokens = tokenize(line);
        
        if (tokens.empty() || tokens[0] == "END" || tokens[0] == "NBFIX") break;
        
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
        skipComments(input);
        auto tokens = tokenize(line);
        
        if (tokens.empty() || tokens[0] == "END") break;
        
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
