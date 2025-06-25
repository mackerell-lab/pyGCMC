// src/io/forcefield/PrmParserOperations.hpp

#pragma once

#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

// Forward declarations - resolved by ModelModule.hpp in implementation files

class PrmParserOperations {
public:
    // Section parsing operations
    static void parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine);
    static void parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff);
    
    // Stream parsing
    static void parseStream(std::istream& input, pygcmc::model::ForceField& ff);
    
    // High-level parsing operations
    static void parse_string(const std::string& content, pygcmc::model::ForceField& ff);
    static void parse_file_to_forcefield(const std::string& filename, pygcmc::model::ForceField& ff);
    static pygcmc::model::ForceField parse_file(const std::string& filename);
    static pygcmc::model::ForceField parse_files(const std::vector<std::string>& filenames);
};

} // namespace io
} // namespace pygcmc