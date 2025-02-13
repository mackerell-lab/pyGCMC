// src/platform/cpu/energy.hpp

#pragma once

#include "../../model/montecarlo.hpp"
#include "../platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
extern const float COULOMB;
extern const float MIN_SAFE_DISTANCE;
extern const float MAX_SAFE_ENERGY;

// Ewald parameters
struct EwaldParams {
    float alpha{1.0f};     // Ewald分离参数 (nm^-1)
    int kmax[3]{6,6,6};    // 倒空间最大波矢
    float tolerance{1e-5f}; // 精度控制
    bool initialized{false};
};

// Global Ewald parameters
extern EwaldParams ewald_params;

/**
 * @brief Calculate nonbonded energies for movement residues only
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between movement residues and all other active residues without distance cutoff.
 * 
 * @param state System state containing residues and force field parameters
 */
void computeMovementEnergy(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for movement residues only with distance cutoff
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between movement residues and all other active residues within the specified cutoff distance.
 * 
 * @param state System state containing residues and force field parameters
 */
void computeMovementEnergyCutoff(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues without distance cutoff.
 * 
 * @param state System state containing residues and force field parameters
 */
void computeSystemEnergy(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with distance cutoff
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues within the specified cutoff distance.
 * 
 * @param state System state containing residues and force field parameters
 */
void computeSystemEnergyCutoff(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with periodic boundary conditions
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues without distance cutoff,
 * applying periodic boundary conditions using the minimum image convention.
 * 
 * @param state System state containing residues and force field parameters
 * @throws std::runtime_error if box dimensions are invalid for PBC calculation
 */
void computeSystemEnergyPBC(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with periodic boundary conditions and cutoff
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues within the specified cutoff distance,
 * applying periodic boundary conditions using the minimum image convention.
 * 
 * @param state System state containing residues and force field parameters
 * @throws std::runtime_error if box dimensions are invalid for PBC calculation
 */
void computeSystemEnergyPBCCutoff(model::MCState& state);

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable);

// Ewald method interfaces
void setEwaldParameters(float alpha, const int kmax[3], float tolerance = 1e-5f);
void computeSystemEnergyEwald(model::MCState& state);
void computeMovementEnergyEwald(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc