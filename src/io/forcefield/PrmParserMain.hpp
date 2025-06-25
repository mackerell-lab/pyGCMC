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
    static void parse_string(const std::string& content, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parse_string(content, ff);
    }
    
    static void parse_file_to_forcefield(const std::string& filename, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parse_file_to_forcefield(filename, ff);
    }
    
    static pygcmc::model::ForceField parse_file(const std::string& filename) {
        return PrmParserOperations::parse_file(filename);
    }
    
    static pygcmc::model::ForceField parse_files(const std::vector<std::string>& filenames) {
        return PrmParserOperations::parse_files(filenames);
    }

    // Instance method for backward compatibility
    void parse(const std::string& filename, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parse_file_to_forcefield(filename, ff);
    }

private:
    // Helper functions (delegate to structures for compatibility)
    void skipComments(std::istream&) {
        // Not used in current implementation, kept for compatibility
    }
    
    std::vector<std::string> tokenize(const std::string& line) {
        return PrmParserStructures::tokenize(line);
    }
    
    // Section detection (delegate to structures)
    bool isAtomsSection(const std::string& line) {
        return PrmParserStructures::isAtomsSection(line);
    }
    bool isBondsSection(const std::string& line) {
        return PrmParserStructures::isBondsSection(line);
    }
    bool isAnglesSection(const std::string& line) {
        return PrmParserStructures::isAnglesSection(line);
    }
    bool isDihedralsSection(const std::string& line) {
        return PrmParserStructures::isDihedralsSection(line);
    }
    bool isImproperSection(const std::string& line) {
        return PrmParserStructures::isImproperSection(line);
    }
    bool isNonbondedSection(const std::string& line) {
        return PrmParserStructures::isNonbondedSection(line);
    }
    bool isNBFixSection(const std::string& line) {
        return PrmParserStructures::isNBFixSection(line);
    }
    
    // Section parsing (delegate to operations)
    void parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parseAtomsSection(input, ff);
    }
    void parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parseBondsSection(input, ff);
    }
    void parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parseAnglesSection(input, ff);
    }
    void parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parseDihedralsSection(input, ff);
    }
    void parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parseImproperSection(input, ff);
    }
    void parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine) {
        PrmParserOperations::parseNonbondedSection(input, ff, firstLine);
    }
    void parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parseNBFixSection(input, ff);
    }
    
    // Stream parsing (delegate to operations)
    void parseStream(std::istream& input, pygcmc::model::ForceField& ff) {
        PrmParserOperations::parseStream(input, ff);
    }
    
    // Line processing helpers (delegate to structures)
    std::string readContinuationLine(std::istream& input, std::string firstLine) {
        return PrmParserStructures::readContinuationLine(input, firstLine);
    }
    bool isCommentLine(const std::string& line) {
        return PrmParserStructures::isCommentLine(line);
    }
    std::string removeComments(const std::string& line) {
        return PrmParserStructures::removeComments(line);
    }
    std::string trim(const std::string& str) {
        return PrmParserStructures::trim(str);
    }

    // Parameter processing helpers (delegate to structures)
    std::pair<std::string, std::string> make_type_pair(const std::string& type1, const std::string& type2) const {
        return PrmParserStructures::make_type_pair(type1, type2);
    }
    std::tuple<std::string, std::string, std::string> make_type_triple(
        const std::string& type1, const std::string& type2, const std::string& type3) const {
        return PrmParserStructures::make_type_triple(type1, type2, type3);
    }
    std::tuple<std::string, std::string, std::string, std::string> make_type_quad(
        const std::string& type1, const std::string& type2, 
        const std::string& type3, const std::string& type4) const {
        return PrmParserStructures::make_type_quad(type1, type2, type3, type4);
    }
};

} // namespace io
} // namespace pygcmc