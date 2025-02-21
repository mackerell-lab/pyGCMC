#include "energyCommon.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
const float COULOMB = 138.935456f;
const float MIN_SAFE_DISTANCE = 0.01f;  // nm (1% of typical sigma)
const float MAX_SAFE_ENERGY = 1e6f;     // kJ/mol

} // namespace cpu
} // namespace platform
} // namespace pygcmc 