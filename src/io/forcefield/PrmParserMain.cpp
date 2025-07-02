// src/io/forcefield/PrmParserMain.cpp

#include "PrmParserMain.hpp"
#include "PrmParserOperations.hpp"
#include "PrmParserSections.hpp"
#include "model/ModelModule.hpp"
#include <fstream>

namespace pygcmc {
namespace io {

// Initialize static debug flag
bool PRMParser::debug_output = false;

void PRMParser::parse_string(const std::string& content, pygcmc::model::ForceField& ff) {
    // Sync debug flags
    PrmParserOperations::getDebugFlag() = debug_output;
    PrmParserSections::getDebugFlag() = debug_output;
    std::istringstream iss(content);
    PrmParserOperations::parseStream(iss, ff);
}

void PRMParser::parse_file_to_forcefield(const std::string& filename, pygcmc::model::ForceField& ff) {
    // Sync debug flags
    PrmParserOperations::getDebugFlag() = debug_output;
    PrmParserSections::getDebugFlag() = debug_output;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open parameter file: " + filename);
    }
    PrmParserOperations::parseStream(file, ff);
}

pygcmc::model::ForceField PRMParser::parse_file(const std::string& filename) {
    // Sync debug flags
    PrmParserOperations::getDebugFlag() = debug_output;
    PrmParserSections::getDebugFlag() = debug_output;
    pygcmc::model::ForceField ff;
    parse_file_to_forcefield(filename, ff);
    return ff;
}

pygcmc::model::ForceField PRMParser::parse_files(const std::vector<std::string>& filenames) {
    pygcmc::model::ForceField ff;
    for (const auto& filename : filenames) {
        parse_file_to_forcefield(filename, ff);
    }
    return ff;
}

// Instance method for backward compatibility
void PRMParser::parse(const std::string& filename, pygcmc::model::ForceField& ff) {
    parse_file_to_forcefield(filename, ff);
}

// Delegate section parsing methods to PrmParserOperations
void PRMParser::parseStream(std::istream& input, pygcmc::model::ForceField& ff) {
    // Sync debug flags
    PrmParserOperations::getDebugFlag() = debug_output;
    PrmParserOperations::parseStream(input, ff);
}

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
    PrmParserSections::parseDihedralsSection(input, ff);
}

void PRMParser::parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserSections::parseImproperSection(input, ff);
}

void PRMParser::parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine) {
    PrmParserSections::parseNonbondedSection(input, ff, firstLine);
}

void PRMParser::parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff) {
    PrmParserSections::parseNBFixSection(input, ff);
}

} // namespace io
} // namespace pygcmc