// src/io/prmParser.hpp

#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include "model/forcefield.hpp"

namespace pygcmc {
namespace io {

class PRMParser {
public:
    // Static debug flag
    static bool debug_output;

    PRMParser() = default;
    ~PRMParser() = default;

    // Static methods for parsing
    static void parse_string(const std::string& content, model::ForceField& ff);
    static void parse_file_to_forcefield(const std::string& filename, model::ForceField& ff);
    static model::ForceField parse_file(const std::string& filename);
    static model::ForceField parse_files(const std::vector<std::string>& filenames);

    // Instance method for backward compatibility
    void parse(const std::string& filename, model::ForceField& ff);

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
    void parseAtomsSection(std::istream& input, model::ForceField& ff);
    void parseBondsSection(std::istream& input, model::ForceField& ff);
    void parseAnglesSection(std::istream& input, model::ForceField& ff);
    void parseDihedralsSection(std::istream& input, model::ForceField& ff);
    void parseImproperSection(std::istream& input, model::ForceField& ff);
    void parseNonbondedSection(std::istream& input, model::ForceField& ff, const std::string& firstLine);
    void parseNBFixSection(std::istream& input, model::ForceField& ff);
    
    // Stream parsing
    void parseStream(std::istream& input, model::ForceField& ff);
    
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

} // namespace io
} // namespace pygcmc
