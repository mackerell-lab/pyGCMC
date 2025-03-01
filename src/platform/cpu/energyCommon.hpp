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

// Debug output control
extern bool energy_debug_output;

// Energy calculation method enumeration
enum class EnergyMethod {
    DIRECT,  // Direct calculation method
    EWALD    // Ewald summation method
};

// CHARMM-style switching function parameters
struct SwitchingFunctionParams {
    bool use_switching{false};  // Whether to use switching function (default: false)
    float r_on{0.8f};           // Inner cutoff radius (nm) where switching starts
    float r_off{1.0f};          // Outer cutoff radius (nm) where potential goes to zero
};

// Global switching function parameters
extern SwitchingFunctionParams switching_params;

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

// Calculate VDW energy (shared by both methods)
inline double calculateVdwEnergy(double r2, double sigma, double eps) {
    double sigma_r2 = (sigma * sigma) / r2;
    double sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;
    double sigma_r12 = sigma_r6 * sigma_r6;
    return 4.0 * eps * (sigma_r12 - sigma_r6);
}

// Calculate CHARMM switching function S(r)
inline float calculateSwitchingFunction(float r) {
    if (!switching_params.use_switching || r <= switching_params.r_on) {
        return 1.0f;  // No switching below r_on
    }
    if (r >= switching_params.r_off) {
        return 0.0f;  // Zero potential beyond r_off
    }
    
    // Calculate CHARMM-style switching function
    // S(r) = [(r_off^2 - r^2)^2 * (r_off^2 + 2r^2 - 3r_on^2)] / (r_off^2 - r_on^2)^3
    float r2 = r * r;
    float ron2 = switching_params.r_on * switching_params.r_on;
    float roff2 = switching_params.r_off * switching_params.r_off;
    
    float numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0f*r2 - 3.0f*ron2);
    float denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2);
    
    return numerator / denominator;
}

// Calculate VDW energy with CHARMM switching function applied
inline double calculateSwitchedVdwEnergy(double r2, double sigma, double eps) {
    float r = std::sqrt(r2);
    
    // Get the basic LJ energy
    double energy = calculateVdwEnergy(r2, sigma, eps);
    
    // Apply switching function if enabled and necessary
    if (switching_params.use_switching) {
        float switch_val = calculateSwitchingFunction(r);
        energy *= switch_val;
    }
    
    return energy;
}

// Configure CHARMM-style switching function
void setSwitchingFunction(bool use_switching, float r_on, float r_off);

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

// Set energy debug output
inline void setEnergyDebugOutput(bool enable) {
    energy_debug_output = enable;
}

// Debug output helper function
inline void logEnergyDebug(const std::string& message) {
    if (energy_debug_output) {
        platform::log(LogLevel::DEBUG, message);
    }
}

// 前置声明所有direct计算中的函数，避免歧义
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

// Declare Direct and Ewald calculation functions (internal use)
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeSystemEnergyEwald(model::MCState& state);
void computeMovementEnergyEwald(model::MCState& state);

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
    // 计算所有residue中的能量（包含vdw能量和实空间静电能量）
    double residue_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            residue_total += residue.energy_vdw + residue.energy_elec;
        }
    }
    
    // 加上EwaldEnergy中的倒空间能量和自能
    return residue_total + state.ewald_energy.reciprocal + state.ewald_energy.self;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 