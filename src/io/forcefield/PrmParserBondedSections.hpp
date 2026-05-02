// src/io/forcefield/PrmParserBondedSections.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

class PrmParserBondedSections {
public:
    // Parse dihedrals section
    static void parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);

    // Parse improper section
    static void parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);

    // Parse NBFIX section
    static void parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);
};

} // namespace io
} // namespace pygcmc
