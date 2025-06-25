// src/io/forcefield/PrmParserMain.hpp

#pragma once

#include "PrmParserStructures.hpp"
#include "PrmParserOperations.hpp"
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace io {

// Forward declarations resolved by ModelModule.hpp

class PRMParser {
public:
    // Static debug flag (delegate to structures for compatibility)
    static bool debug_output;

    PRMParser() = default;
    ~PRMParser() = default;

    // Static methods for parsing (delegate to operations)
    static void parse_string(const std::string& content, pygcmc::model::ForceField& ff);
    static void parse_file_to_forcefield(const std::string& filename, pygcmc::model::ForceField& ff);
    static pygcmc::model::ForceField parse_file(const std::string& filename);
    static pygcmc::model::ForceField parse_files(const std::vector<std::string>& filenames);

    // Instance method for backward compatibility
    void parse(const std::string& filename, pygcmc::model::ForceField& ff);

private:
    // Helper functions (delegate to structures for compatibility)
    void skipComments(std::istream&);
    std::vector<std::string> tokenize(const std::string& line);
    
    // Section detection (delegate to structures)
    bool isAtomsSection(const std::string& line);
    bool isBondsSection(const std::string& line);
    bool isAnglesSection(const std::string& line);
    bool isDihedralsSection(const std::string& line);
    bool isImproperSection(const std::string& line);
    bool isNonbondedSection(const std::string& line);
    bool isNBFixSection(const std::string& line);
    
    // Section parsing (delegate to operations)
    void parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff);
    void parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff);
    void parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff);
    void parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff);
    void parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff);
    void parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine);
    void parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff);
    
    // Stream parsing (delegate to operations)
    void parseStream(std::istream& input, pygcmc::model::ForceField& ff);
    
    // Line processing helpers (delegate to structures)
    std::string readContinuationLine(std::istream& input, std::string firstLine);
    bool isCommentLine(const std::string& line);
    std::string removeComments(const std::string& line);
    std::string trim(const std::string& str);

    // Parameter processing helpers (delegate to structures)
    std::pair<std::string, std::string> make_type_pair(const std::string& type1, const std::string& type2) const;
    std::tuple<std::string, std::string, std::string> make_type_triple(
        const std::string& type1, const std::string& type2, const std::string& type3) const;
    std::tuple<std::string, std::string, std::string, std::string> make_type_quad(
        const std::string& type1, const std::string& type2, 
        const std::string& type3, const std::string& type4) const;
};

} // namespace io
} // namespace pygcmc