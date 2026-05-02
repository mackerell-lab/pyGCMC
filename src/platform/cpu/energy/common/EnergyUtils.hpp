#pragma once

#include "EnergyConstants.hpp"
#include "model/ModelModule.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Debug output control
extern bool energy_debug_output;

// Function to get energy debug output status (using platform's debug_mode)
inline bool getEnergyDebugOutput() {
    return platform::is_debug_mode();
}

// Common utility functions
inline float capEnergy(float energy) {
    return std::min(std::max(energy, -MAX_SAFE_ENERGY), MAX_SAFE_ENERGY);
}

inline float checkDistance(float r2) {
    const float min_r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;
    return (r2 < min_r2) ? min_r2 : r2;
}

// Apply periodic boundary conditions
inline void applyPBC(float& dx, float& dy, float& dz, const float box[3]) {
    if(dx > box[0]/2) dx -= box[0];
    else if(dx < -box[0]/2) dx += box[0];
    if(dy > box[1]/2) dy -= box[1];
    else if(dy < -box[1]/2) dy += box[1];
    if(dz > box[2]/2) dz -= box[2];
    else if(dz < -box[2]/2) dz += box[2];
}

// Verify PBC box
inline void validateBox(const float box[3], float cutoff = 0.0f) {
    if (box[0] <= 0.0f || box[1] <= 0.0f || box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }

    if (cutoff > 0.0f) {
        float minBoxSize = std::min(box[0], std::min(box[1], box[2]));
        if (cutoff >= 0.5f * minBoxSize) {
            platform::log(LogLevel::WARNING,
                "Warning: Cutoff distance (", cutoff,
                " nm) is larger than half the smallest box dimension (",
                minBoxSize/2, " nm). This may affect minimum image convention.");
        }
    }
}

// Verify system neutrality
inline void checkSystemNeutrality(const model::MCState& state) {
    float totalCharge = 0.0f;
    for(const auto& atom : state.atoms) {
        totalCharge += atom.charge;
    }
    if(std::abs(totalCharge) > 1e-6f) {
        throw std::runtime_error("Energy calculation requires neutral system");
    }
}

// Replace setEnergyDebugOutput function
inline void setEnergyDebugOutput(bool enable) {
    energy_debug_output = enable;
    // Also update the platform's debug mode for consistency
    platform::set_debug_mode(enable);
}

inline void logEnergyDebug(const std::string& message) {
    if (getEnergyDebugOutput()) {
        platform::log(LogLevel::DEBUG, message);
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
