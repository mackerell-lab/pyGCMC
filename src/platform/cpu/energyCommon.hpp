// src/platform/cpu/energyCommon.hpp

#pragma once

#include <cmath>
#include "model/montecarlo.hpp"
#include "platform/platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
extern const float COULOMB;      // Coulomb constant in GROMACS MD units [kJ·nm/mol/e²]
extern const float MIN_SAFE_DISTANCE;  // Minimum allowed distance (1% of sigma)
extern const float MAX_SAFE_ENERGY;    // Maximum allowed energy per interaction

// Constants for Ewald and PME calculations
static const int NUM_TABLE_POINTS = 20000;  // High precision table size for approximations
static const double TWO_OVER_SQRT_PI = 2.0/std::sqrt(M_PI);  // Constant for Ewald calculations

// Debug output control
extern bool energy_debug_output;

// Function to get energy debug output status (using platform's debug_mode)
inline bool getEnergyDebugOutput() {
    return platform::is_debug_mode();
}

// Energy calculation method enumeration
enum class EnergyMethod {
    DIRECT,  // Direct calculation method
    EWALD,   // Ewald summation method
    PME      // Particle Mesh Ewald method
};

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

// NOTE: calculateSwitchingFunction 函数已移至 energyLJ.hpp 和 energyLJ.cpp 中
// 请使用 platform::cpu::calculateSwitchingFunction 函数

// Configure CHARMM-style switching function
void setSwitchingFunction(model::MCState& state, bool use_switching, float r_on, float r_off);

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

// Forward declarations for all direct calculation functions to avoid ambiguity
void computeMovementEnergy(model::MCState& state);
void computeMovementEnergyCutoff(model::MCState& state);
void computeSystemEnergy(model::MCState& state);
void computeSystemEnergyCutoff(model::MCState& state);
void computeSystemEnergyPBC(model::MCState& state);
void computeSystemEnergyPBCCutoff(model::MCState& state);
void computeSystemVdwEnergyCutoff(model::MCState& state);

// Unified energy calculation interface
void computeSystemEnergy(model::MCState& state, 
                         EnergyMethod method = EnergyMethod::DIRECT,
                         bool use_cutoff = false, 
                         bool use_pbc = false);

void computeMovementEnergy(model::MCState& state, 
                          EnergyMethod method = EnergyMethod::DIRECT,
                          bool use_cutoff = false, 
                          bool use_pbc = false);

// Declare Direct, Ewald and PME calculation functions (internal use)
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeSystemEnergyEwald(model::MCState& state);
void computeMovementEnergyEwald(model::MCState& state);
void computeSystemEnergyPME(model::MCState& state);
void computeMovementEnergyPME(model::MCState& state);

// Get total energy function
inline double getTotalEnergy(const model::MCState& state) {
    double total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            total += residue.energy_vdw + residue.energy_elec;
        }
    }
    return total;
}

inline double getEwaldTotalEnergy(const model::MCState& state) {
    // Calculate energy from all residues (including vdw energy and real-space electrostatic energy)
    double residue_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            residue_total += residue.energy_vdw + residue.energy_elec;
        }
    }
    
    // Add reciprocal-space energy and self-energy from EwaldEnergy
    return residue_total + state.ewald_energy.reciprocal + state.ewald_energy.self;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 