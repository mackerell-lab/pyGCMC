// src/io/forcefield/PrmParserDrudeScan.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

class PrmParserDrudeScan {
public:
    // Pre-scan the input stream for Drude-related parameters
    static void prescanForDrudeParameters(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output);
};

} // namespace io
} // namespace pygcmc