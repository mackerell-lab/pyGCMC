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
    // Parse nonbonded section (the most complex section)
    static void parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine);
    
    // Debug flag access
    static bool& getDebugFlag();
};

} // namespace io
} // namespace pygcmc