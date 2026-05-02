// src/io/forcefield/PrmParserOperations.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

class PrmParserOperations {
public:
    // Stream parsing - the main parsing loop
    static void parseStream(std::istream& input, pygcmc::model::ForceField& ff);

    // Debug flag access
    static bool& getDebugFlag();
};

} // namespace io
} // namespace pygcmc
