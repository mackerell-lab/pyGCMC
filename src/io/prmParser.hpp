#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include "../model/forcefield.hpp"

namespace pygcmc {

class PrmParser {
public:
    PrmParser() = default;
    ~PrmParser() = default;

    // Static methods for parsing
    static void parse_string(const std::string& content, ForceField& ff);
    static void parse_file(const std::string& filename, ForceField& ff);

    // Instance method for backward compatibility
    void parse(const std::string& filename, ForceField& ff);

private:
    // Helper functions
    void skipComments(std::istream& input);
    std::vector<std::string> tokenize(const std::string& line);
    bool isNonbondedSection(const std::string& line);
    bool isNBFixSection(const std::string& line);
    void parseNonbondedSection(std::istream& input, ForceField& ff, const std::string& firstLine);
    void parseNBFixSection(std::istream& input, ForceField& ff);
    void parseStream(std::istream& input, ForceField& ff);
};

} // namespace pygcmc
