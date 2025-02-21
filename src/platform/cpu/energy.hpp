#pragma once

// Include both energy calculation implementations
#include "energyDirect.hpp"
#include "energyEwald.hpp"

// This header serves as a unified interface for all energy calculation methods
// All declarations are now in the respective headers:
// - energyDirect.hpp: Direct calculation methods and constants
// - energyEwald.hpp: Ewald summation methods and parameters 