// src/io/prmParser.hpp

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
    
    // Section detection
    bool isAtomsSection(const std::string& line);
    bool isBondsSection(const std::string& line);
    bool isAnglesSection(const std::string& line);
    bool isDihedralsSection(const std::string& line);
    bool isImproperSection(const std::string& line);
    bool isNonbondedSection(const std::string& line);
    bool isNBFixSection(const std::string& line);
    
    // Section parsing
    void parseAtomsSection(std::istream& input, ForceField& ff);
    void parseBondsSection(std::istream& input, ForceField& ff);
    void parseAnglesSection(std::istream& input, ForceField& ff);
    void parseDihedralsSection(std::istream& input, ForceField& ff);
    void parseImproperSection(std::istream& input, ForceField& ff);
    void parseNonbondedSection(std::istream& input, ForceField& ff, const std::string& firstLine);
    void parseNBFixSection(std::istream& input, ForceField& ff);
    
    // Stream parsing
    void parseStream(std::istream& input, ForceField& ff);
    
    // Line processing helpers
    std::string readContinuationLine(std::istream& input, std::string firstLine);
    bool isCommentLine(const std::string& line);
    std::string removeComments(const std::string& line);
    std::string trim(const std::string& str);

    // Parameter processing helpers
    std::pair<std::string, std::string> make_type_pair(const std::string& type1, const std::string& type2) const;
    std::tuple<std::string, std::string, std::string> make_type_triple(
        const std::string& type1, const std::string& type2, const std::string& type3) const;
    std::tuple<std::string, std::string, std::string, std::string> make_type_quad(
        const std::string& type1, const std::string& type2, 
        const std::string& type3, const std::string& type4) const;
};

} // namespace pygcmc
