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
std::mt19937 RotationUtils::gen_(std::random_device{}());
std::uniform_real_distribution<> RotationUtils::dis_(0.0, 1.0);

} // namespace utils
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc