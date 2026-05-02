// src/io/forcefield/PrmParserDrudeSections.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

class PrmParserDrudeSections {
public:
    // Parse ALPHA/THOLE section
    static void parseAlphaTHoleSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);

    // Parse LONEPAIR section
    static void parseLonePairSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);

    // Parse ANISOTROPY section
    static void parseAnisotropySection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);
};

} // namespace io
} // namespace pygcmc
