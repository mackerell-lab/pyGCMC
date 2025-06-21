// src/model/montecarlo.hpp

#pragma once

// Include the new refactored montecarlo module
#include "montecarlo/MCMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using MCState = montecarlo::MCState;
using MCResidue = montecarlo::MCResidue;
using MCAtom = montecarlo::MCAtom;

} // namespace model
} // namespace pygcmc
