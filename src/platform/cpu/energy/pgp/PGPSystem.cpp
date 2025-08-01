#include "PGPSystem.hpp"
#include "PGPGlobal.hpp"
#include "PGPCore.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"
#include "PGPInterpolation.hpp"
#include "PGPPrecompute.hpp"
#include "../lj/LJMain.hpp"
#include "../common/EnergyDirectCore.hpp"
#include "platform/platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Use PGP method to calculate system energy
 */
void computeSystemEnergyPGPImpl(model::MCState& state) {
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing total system energy using PGP method");
    }
    
    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate real space part
    computeRealSpacePGPImpl(state, false, true);
    
    // 3. Calculate self energy correction
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, false);
    
    // 4. Calculate LJ interactions using direct cutoff method
    computeSystemVdwEnergyCutoff(state);
    
    // Multiply real space energy by COULOMB constant
    state.ewald_energy.real_space *= COULOMB;
    
    // Calculate total LJ energy
    double vdw_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            vdw_total += residue.energy_vdw;
        }
    }
    
    // Calculate total energy
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP system energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

/**
 * @brief Use PGP method to calculate energy of moving residues
 */
void computeMovementEnergyPGPImpl(model::MCState& state) {
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing movement residue energy using PGP method");
    }
    
    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate real space part (only for moving residues)
    computeRealSpacePGPImpl(state, true, true);
    
    // 3. Calculate self energy correction (only for moving residues)
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, true);
    
    // 4. Calculate LJ interactions using direct cutoff method
    computeSystemVdwEnergyCutoff(state);
    
    // Multiply real space energy by COULOMB constant
    state.ewald_energy.real_space *= COULOMB;
    
    // Only accumulate LJ energy for moving residues
    double vdw_total = 0.0;
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                vdw_total += state.residues[i].energy_vdw;
            }
        }
    }
    
    // Calculate total energy
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP movement energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

/**
 * @brief Public interface wrapper for system energy calculation
 */
void computeSystemEnergyPGP(model::MCState& state) {
    computeSystemEnergyPGPImpl(state);
}

/**
 * @brief Public interface wrapper for movement energy calculation
 */
void computeMovementEnergyPGP(model::MCState& state) {
    computeMovementEnergyPGPImpl(state);
}



} // namespace cpu
} // namespace platform
} // namespace pygcmc 