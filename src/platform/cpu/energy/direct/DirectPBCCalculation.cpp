#include "DirectPBCCalculation.hpp"
#include "DirectCore.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <stdexcept>
#include <sstream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

void computeSystemEnergyPBC(model::MCState& state) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting PBC nonbonded energy calculation (no cutoff) ===";
        ss << "\nBox dimensions: " << state.info.box[0] << " x " 
           << state.info.box[1] << " x " << state.info.box[2] << " nm";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Validate box dimensions before proceeding
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    // Calculate with PBC enabled but no cutoff
    computeNonbondedEnergy(state, false, false, true);
}

void computeSystemEnergyPBCCutoff(model::MCState& state) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting PBC nonbonded energy calculation (with cutoff) ===";
        ss << "\nBox dimensions: " << state.info.box[0] << " x " 
           << state.info.box[1] << " x " << state.info.box[2] << " nm";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Validate box dimensions before proceeding
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    // Calculate with both cutoff and PBC enabled
    computeNonbondedEnergy(state, true, false, true);
}

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 