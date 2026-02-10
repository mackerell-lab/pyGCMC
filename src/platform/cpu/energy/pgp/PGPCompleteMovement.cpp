// PGPCompleteMovement.cpp
// Correct implementation of PGP Complete for movement energy calculation

#include "PGPComplete.hpp"
#include "PGPGlobal.hpp"
#include "PGPCore.hpp"
#include "PGPInterpolation.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Helper function to check if an atom belongs to movement residues
static bool isMovementAtom(int atomIndex, const model::MCState& state) {
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            const auto& res = state.residues[i];
            if (res.active && atomIndex >= res.atomStart && 
                atomIndex < res.atomStart + res.atomCount) {
                return true;
            }
        }
    }
    return false;
}

// Calculate real space electrostatics for movement atoms
static void calculateMovementRealSpace(model::MCState& state) {
    const auto& atoms = state.atoms;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const double alpha = getPGPParams().alpha;
    
    // Reset real space energy
    state.ewald_energy.real_space = 0.0;
    
    // For each movement atom, calculate interactions with ALL atoms
    for (int i = 0; i < state.activeAtomCount; ++i) {
        if (!isMovementAtom(i, state)) continue;
        
        const double qi = atoms[i].charge;
        if (qi == 0.0) continue;
        
        const double xi = atoms[i].x;
        const double yi = atoms[i].y;
        const double zi = atoms[i].z;
        
        // Interact with all other atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            if (i == j) continue;  // Skip self
            
            const double qj = atoms[j].charge;
            if (qj == 0.0) continue;
            
            // Calculate distance with PBC
            double dx = xi - atoms[j].x;
            double dy = yi - atoms[j].y;
            double dz = zi - atoms[j].z;
            
            // Apply minimum image convention
            float dx_f = static_cast<float>(dx);
            float dy_f = static_cast<float>(dy);
            float dz_f = static_cast<float>(dz);
            applyPBC(dx_f, dy_f, dz_f, state.info.box);
            dx = dx_f; dy = dy_f; dz = dz_f;
            
            const double r2 = dx*dx + dy*dy + dz*dz;
            
            // Apply cutoff
            if (r2 > cutoff2) continue;
            
            // Skip extremely close atoms
            if (r2 < 1e-12) continue;
            
            const double r = std::sqrt(r2);
            
            // Calculate erfc(alpha*r)/r
            const double alphar = alpha * r;
            const double erfc_val = std::erfc(alphar);
            
            // Real space electrostatic energy
            // Note: We count each movement atom's interaction once
            const double energy = qi * qj * erfc_val / r;
            
            state.ewald_energy.real_space += energy;
        }
    }
}

// Calculate LJ energy for movement atoms
static void calculateMovementLJ(model::MCState& state) {
    const auto& atoms = state.atoms;
    const auto& forcefield = state.forcefield;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    
    // Reset VDW energies for movement residues
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                state.residues[i].energy_vdw = 0.0f;
            }
        }
    }
    
    // For each movement atom, calculate LJ with ALL atoms
    for (int i = 0; i < state.activeAtomCount; ++i) {
        if (!isMovementAtom(i, state)) continue;
        
        const int type_i = atoms[i].type;
        const double xi = atoms[i].x;
        const double yi = atoms[i].y;
        const double zi = atoms[i].z;
        
        // Find which residue this atom belongs to
        int res_i = -1;
        for (int r = 0; r < state.activeResidueCount; ++r) {
            const auto& res = state.residues[r];
            if (res.active && i >= res.atomStart && i < res.atomStart + res.atomCount) {
                res_i = r;
                break;
            }
        }
        
        if (res_i == -1) continue;
        
        // Interact with all atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            if (i == j) continue;  // Skip self
            
            const int type_j = atoms[j].type;
            
            // Get LJ parameters
            int param_idx = type_i * forcefield.numTotalTypes + type_j;
            const double sigma = forcefield.ljSigma[param_idx];
            const double epsilon = forcefield.ljEps[param_idx];
            
            if (epsilon == 0.0 || sigma == 0.0) continue;
            
            // Calculate distance with PBC
            double dx = xi - atoms[j].x;
            double dy = yi - atoms[j].y;
            double dz = zi - atoms[j].z;
            
            // Apply minimum image convention
            float dx_f = static_cast<float>(dx);
            float dy_f = static_cast<float>(dy);
            float dz_f = static_cast<float>(dz);
            applyPBC(dx_f, dy_f, dz_f, state.info.box);
            dx = dx_f; dy = dy_f; dz = dz_f;
            
            const double r2 = dx*dx + dy*dy + dz*dz;
            
            // Apply cutoff
            if (r2 > cutoff2) continue;
            
            // Skip extremely close atoms
            if (r2 < 1e-12) continue;
            
            // Calculate LJ energy
            const double sigma2 = sigma * sigma;
            const double sigma6 = sigma2 * sigma2 * sigma2;
            const double sigma12 = sigma6 * sigma6;
            const double r6 = r2 * r2 * r2;
            const double r12 = r6 * r6;
            
            const double lj_energy = 4.0 * epsilon * (sigma12/r12 - sigma6/r6);
            
            // Add to residue energy
            state.residues[res_i].energy_vdw += lj_energy;
        }
    }
}

// Calculate self energy for movement atoms only
static double calculateMovementSelfEnergy(const model::MCState& state) {
    const double alpha = getPGPParams().alpha;
    const double self_factor = -alpha / std::sqrt(M_PI) * COULOMB;
    
    double self_energy = 0.0;
    
    for (int i = 0; i < state.activeAtomCount; ++i) {
        if (!isMovementAtom(i, state)) continue;
        
        const double q = state.atoms[i].charge;
        self_energy += q * q * self_factor;
    }
    
    return self_energy;
}

// Legacy implementation kept for reference only.
// The active corrected implementation is in PGPCompleteCorrected.cpp.
void computeMovementEnergyPGPCompleteCorrectLegacy(model::MCState& state) {
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    platform::log(LogLevel::DEBUG, "Computing movement energy using corrected PGP Complete method");
    
    // 1. Calculate grid potential interpolation (reciprocal space contribution)
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate real space for movement atoms
    calculateMovementRealSpace(state);
    
    // 3. Calculate self energy for movement atoms
    state.ewald_energy.self = calculateMovementSelfEnergy(state);
    
    // 4. Calculate LJ for movement atoms
    calculateMovementLJ(state);
    
    // Apply COULOMB constant to real space
    state.ewald_energy.real_space *= COULOMB;
    
    // Calculate total VDW from movement residues
    double vdw_total = 0.0;
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                vdw_total += state.residues[i].energy_vdw;
            }
        }
    }
    
    // Set reciprocal and total
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    platform::log(LogLevel::INFO, "PGP Complete Movement (corrected):");
    platform::log(LogLevel::INFO, "  Grid (reciprocal): ", grid_energy);
    platform::log(LogLevel::INFO, "  Real space: ", state.ewald_energy.real_space);
    platform::log(LogLevel::INFO, "  Self: ", state.ewald_energy.self);
    platform::log(LogLevel::INFO, "  VDW: ", vdw_total);
    platform::log(LogLevel::INFO, "  Total: ", state.ewald_energy.total);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
