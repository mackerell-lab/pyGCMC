// src/io/forcefield/PrmParserOperations.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

class PrmParserOperations {
public:
    // Basic section parsing operations
    static void parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff);
    static void parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff);
    
    // Stream parsing - the main parsing loop
    static void parseStream(std::istream& input, pygcmc::model::ForceField& ff);
    
    // Debug flag access
    static bool& getDebugFlag();
};

} // namespace io
} // namespace pygcmc