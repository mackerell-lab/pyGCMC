// src/io/forcefield/PrmParserSections.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

// Forward declaration
class PrmParserOperations;

class PrmParserSections {
public:
    // Simple section parsing operations
    static void parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff);
    
    // Complex section parsing operations  
    static void parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine);
    
    // Drude-specific section parsing operations
    static void parseAlphaTHoleSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseLonePairSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseAnisotropySection(std::istream& input, pygcmc::model::ForceField& ff);
    
    // Debug flag access
    static bool& getDebugFlag();
};

} // namespace io
} // namespace pygcmc