// src/io/forcefield/PrmParserSectionHandlers.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

class PrmParserSectionHandlers {
public:
    // Section parsing operations
    static void parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);
    static void parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);
    static void parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);
};

} // namespace io
} // namespace pygcmc