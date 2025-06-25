// src/io/forcefield/PrmParserMain.cpp

#include "PrmParserMain.hpp"

namespace pygcmc {
namespace io {

// Initialize static debug flag
bool PRMParser::debug_output = false;

// Static methods for parsing (delegate to operations)
void PRMParser::parse_string(const std::string& content, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parse_string(content, ff);
}

void PRMParser::parse_file_to_forcefield(const std::string& filename, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parse_file_to_forcefield(filename, ff);
}

pygcmc::model::ForceField PRMParser::parse_file(const std::string& filename) {
    return PrmParserOperations::parse_file(filename);
}

pygcmc::model::ForceField PRMParser::parse_files(const std::vector<std::string>& filenames) {
    return PrmParserOperations::parse_files(filenames);
}

// Instance method for backward compatibility
void PRMParser::parse(const std::string& filename, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parse_file_to_forcefield(filename, ff);
}

// Helper functions (delegate to structures for compatibility)
void PRMParser::skipComments(std::istream&) {
    // Not used in current implementation, kept for compatibility
}

std::vector<std::string> PRMParser::tokenize(const std::string& line) {
    return PrmParserStructures::tokenize(line);
}

// Section detection (delegate to structures)
bool PRMParser::isAtomsSection(const std::string& line) {
    return PrmParserStructures::isAtomsSection(line);
}

bool PRMParser::isBondsSection(const std::string& line) {
    return PrmParserStructures::isBondsSection(line);
}

bool PRMParser::isAnglesSection(const std::string& line) {
    return PrmParserStructures::isAnglesSection(line);
}

bool PRMParser::isDihedralsSection(const std::string& line) {
    return PrmParserStructures::isDihedralsSection(line);
}

bool PRMParser::isImproperSection(const std::string& line) {
    return PrmParserStructures::isImproperSection(line);
}

bool PRMParser::isNonbondedSection(const std::string& line) {
    return PrmParserStructures::isNonbondedSection(line);
}

bool PRMParser::isNBFixSection(const std::string& line) {
    return PrmParserStructures::isNBFixSection(line);
}

// Section parsing (delegate to operations)
void PRMParser::parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parseAtomsSection(input, ff);
}

void PRMParser::parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parseBondsSection(input, ff);
}

void PRMParser::parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parseAnglesSection(input, ff);
}

void PRMParser::parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parseDihedralsSection(input, ff);
}

void PRMParser::parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parseImproperSection(input, ff);
}

void PRMParser::parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine) {
    PrmParserOperations::parseNonbondedSection(input, ff, firstLine);
}

void PRMParser::parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parseNBFixSection(input, ff);
}

// Stream parsing (delegate to operations)
void PRMParser::parseStream(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserOperations::parseStream(input, ff);
}

// Line processing helpers (delegate to structures)
std::string PRMParser::readContinuationLine(std::istream& input, std::string firstLine) {
    return PrmParserStructures::readContinuationLine(input, firstLine);
}

bool PRMParser::isCommentLine(const std::string& line) {
    return PrmParserStructures::isCommentLine(line);
}

std::string PRMParser::removeComments(const std::string& line) {
    return PrmParserStructures::removeComments(line);
}

std::string PRMParser::trim(const std::string& str) {
    return PrmParserStructures::trim(str);
}

// Parameter processing helpers (delegate to structures)
std::pair<std::string, std::string> PRMParser::make_type_pair(const std::string& type1, const std::string& type2) const {
    return PrmParserStructures::make_type_pair(type1, type2);
}

std::tuple<std::string, std::string, std::string> PRMParser::make_type_triple(
    const std::string& type1, const std::string& type2, const std::string& type3) const {
    return PrmParserStructures::make_type_triple(type1, type2, type3);
}

std::tuple<std::string, std::string, std::string, std::string> PRMParser::make_type_quad(
    const std::string& type1, const std::string& type2, 
    const std::string& type3, const std::string& type4) const {
    return PrmParserStructures::make_type_quad(type1, type2, type3, type4);
}

} // namespace io
} // namespace pygcmc