/**
 * Static definitions for MovementUtils classes
 */

#include "MovementUtils.hpp"
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace utils {

// Static member definitions for RotationUtils
// Initialize with deterministic seed (0), will be set via setSeed()
std::mt19937 RotationUtils::gen_(0u);
std::uniform_real_distribution<> RotationUtils::dis_(0.0, 1.0);

} // namespace utils
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
