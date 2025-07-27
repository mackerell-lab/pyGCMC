// PGPCompleteCorrected.cpp
// Correct implementation of PGP Complete based on PGP principles

#include "PGPComplete.hpp"
#include "PGPCore.hpp"
#include "PGPInterpolation.hpp"
#include "PGPSelf.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <cstdio>

namespace pygcmc {
namespace platform {
namespace cpu {

// Helper function to check if a residue is in movement list
static bool isMovementResidue(int res_idx, const model::MCState& state) {
    // If no movement residues specified, treat all non-fixed residues as movement
    if (state.movementResidues.empty()) {
        return !state.residues[res_idx].fixed;
    }
    
    for (const auto& movementInfo : state.movementResidues) {
        if (res_idx >= movementInfo.startIndex &&
            res_idx < movementInfo.startIndex + movementInfo.activeCount) {
            return true;
        }
    }
    return false;
}

// Calculate real space electrostatics following PGP principles
static void calculateRealSpacePGPComplete(model::MCState& state, bool movement_only) {
    const auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float* box = state.info.box;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const double alpha = pgp_params.alpha;
    
    // Reset real space energy
    state.ewald_energy.real_space = 0.0;
    
    // Reset electrostatic energies in residues
    for (auto& res : residues) {
        if (res.active) {
            res.energy_elec = 0.0f;
        }
    }
    
    double real_space_total = 0.0;
    
    // Loop over residue pairs
    for (int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if (!residues[r1].active) continue;
        
        // For movement_only mode, skip if r1 is not movement
        if (movement_only && !isMovementResidue(r1, state)) continue;
        
        // Start from r1 to include intra-residue interactions
        for (int r2 = r1; r2 < state.activeResidueCount; r2++) {
            if (!residues[r2].active) continue;
            
            // For movement_only mode, at least one residue must be movement
            if (movement_only && !isMovementResidue(r1, state) && !isMovementResidue(r2, state)) {
                continue;
            }
            
            // CRITICAL: If both residues are fixed, skip
            // Their interaction is already in the precomputed grid
            if (residues[r1].fixed && residues[r2].fixed) continue;
            
            // Loop over atom pairs
            for (int i = residues[r1].atomStart; 
                 i < residues[r1].atomStart + residues[r1].atomCount; i++) {
                
                if (i >= state.activeAtomCount) continue;
                
                const double qi = atoms[i].charge;
                if (std::abs(qi) < 1e-6) continue;
                
                int j_start = (r1 == r2) ? i + 1 : residues[r2].atomStart;
                int j_end = residues[r2].atomStart + residues[r2].atomCount;
                
                for (int j = j_start; j < j_end; j++) {
                    if (j >= state.activeAtomCount) continue;
                    
                    const double qj = atoms[j].charge;
                    if (std::abs(qj) < 1e-6) continue;
                    
                    // Calculate distance with PBC
                    double dx = atoms[i].x - atoms[j].x;
                    double dy = atoms[i].y - atoms[j].y;
                    double dz = atoms[i].z - atoms[j].z;
                    
                    // Apply minimum image convention
                    float dx_f = static_cast<float>(dx);
                    float dy_f = static_cast<float>(dy);
                    float dz_f = static_cast<float>(dz);
                    applyPBC(dx_f, dy_f, dz_f, box);
                    dx = dx_f; dy = dy_f; dz = dz_f;
                    
                    const double r2_dist = dx*dx + dy*dy + dz*dz;
                    
                    // Apply cutoff
                    if (r2_dist > cutoff2) continue;
                    
                    // Skip extremely close atoms
                    if (r2_dist < 1e-12) continue;
                    
                    const double r = std::sqrt(r2_dist);
                    
                    // Calculate erfc(alpha*r)/r
                    const double alphar = alpha * r;
                    const double erfc_val = std::erfc(alphar);
                    const double pair_energy = qi * qj * erfc_val / r;
                    
                    // Accumulate energy
                    real_space_total += pair_energy;
                    
                    // Store in residues
                    if (r1 == r2) {
                        // Intra-residue: all energy to this residue
                        residues[r1].energy_elec += pair_energy;
                    } else {
                        // Inter-residue: split energy
                        residues[r1].energy_elec += pair_energy * 0.5;
                        residues[r2].energy_elec += pair_energy * 0.5;
                    }
                }
            }
        }
    }
    
    state.ewald_energy.real_space = real_space_total;
}

// Calculate LJ energy for PGP Complete
static void calculateLJPGPComplete(model::MCState& state, bool movement_only) {
    const auto& atoms = state.atoms;
    const auto& forcefield = state.forcefield;
    auto& residues = state.residues;
    const float* box = state.info.box;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    
    platform::log(LogLevel::DEBUG, "calculateLJPGPComplete: movement_only=", movement_only);
    platform::log(LogLevel::DEBUG, "Active residue count: ", state.activeResidueCount);
    
    // Debug force field
    fprintf(stderr, "DEBUG: calculateLJPGPComplete called\n");
    fprintf(stderr, "  movement_only = %d\n", movement_only);
    fprintf(stderr, "  activeResidueCount = %d\n", state.activeResidueCount);
    fprintf(stderr, "  activeAtomCount = %d\n", state.activeAtomCount);
    fprintf(stderr, "  numTotalTypes = %d\n", forcefield.numTotalTypes);
    fprintf(stderr, "  ljSigma size = %zu\n", forcefield.ljSigma.size());
    fprintf(stderr, "  ljEps size = %zu\n", forcefield.ljEps.size());
    fprintf(stderr, "  cutoff = %f, cutoff2 = %f\n", state.info.cutoff, cutoff2);
    
    // Reset VDW energies
    for (auto& res : residues) {
        if (res.active) {
            res.energy_vdw = 0.0f;
        }
    }
    
    // Loop over residue pairs
    for (int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if (!residues[r1].active) continue;
        
        bool r1_is_movement = isMovementResidue(r1, state);
        printf("  Residue %d: active=%d, fixed=%d, is_movement=%d\n", 
               r1, residues[r1].active, residues[r1].fixed, r1_is_movement);
        
        // For movement_only mode, skip if r1 is not movement
        if (movement_only && !r1_is_movement) {
            platform::log(LogLevel::DEBUG, "Skipping residue ", r1, " - not movement");
            printf("  Skipping residue %d (not movement)\n", r1);
            continue;
        }
        
        platform::log(LogLevel::DEBUG, "Processing residue r1=", r1);
        printf("  Processing residue r1=%d\n", r1);
        
        // Start from r1 to include intra-residue interactions
        for (int r2 = r1; r2 < state.activeResidueCount; r2++) {
            if (!residues[r2].active) continue;
            
            // For movement_only mode, at least one residue must be movement
            if (movement_only && !isMovementResidue(r1, state) && !isMovementResidue(r2, state)) {
                continue;
            }
            
            platform::log(LogLevel::DEBUG, "Processing residue pair r1=", r1, " r2=", r2);
            printf("    Processing residue pair r1=%d, r2=%d\n", r1, r2);
            
            // Loop over atom pairs
            for (int i = residues[r1].atomStart; 
                 i < residues[r1].atomStart + residues[r1].atomCount; i++) {
                
                if (i >= state.activeAtomCount) continue;
                const int type_i = atoms[i].type;
                
                int j_start = (r1 == r2) ? i + 1 : residues[r2].atomStart;
                int j_end = residues[r2].atomStart + residues[r2].atomCount;
                
                for (int j = j_start; j < j_end; j++) {
                    if (j >= state.activeAtomCount) continue;
                    const int type_j = atoms[j].type;
                    
                    // Get LJ parameters
                    int param_idx = type_i * forcefield.numTotalTypes + type_j;
                    
                    printf("      Atom pair i=%d (type %d), j=%d (type %d), param_idx=%d\n", 
                           i, type_i, j, type_j, param_idx);
                    
                    if (param_idx >= static_cast<int>(forcefield.ljSigma.size())) {
                        printf("      ERROR: param_idx %d >= ljSigma size %zu\n", 
                               param_idx, forcefield.ljSigma.size());
                        continue;
                    }
                    
                    const double sigma = forcefield.ljSigma[param_idx];
                    const double epsilon = forcefield.ljEps[param_idx];
                    
                    printf("      LJ params: sigma=%f, epsilon=%f\n", sigma, epsilon);
                    
                    if (epsilon == 0.0 || sigma == 0.0) {
                        printf("      Skipping: epsilon or sigma is 0\n");
                        continue;
                    }
                    
                    // Calculate distance with PBC
                    double dx = atoms[i].x - atoms[j].x;
                    double dy = atoms[i].y - atoms[j].y;
                    double dz = atoms[i].z - atoms[j].z;
                    
                    // Apply minimum image convention
                    float dx_f = static_cast<float>(dx);
                    float dy_f = static_cast<float>(dy);
                    float dz_f = static_cast<float>(dz);
                    applyPBC(dx_f, dy_f, dz_f, box);
                    dx = dx_f; dy = dy_f; dz = dz_f;
                    
                    const double r2_dist = dx*dx + dy*dy + dz*dz;
                    
                    // Apply cutoff
                    if (r2_dist > cutoff2) continue;
                    
                    // Skip extremely close atoms
                    if (r2_dist < 1e-12) continue;
                    
                    // Calculate LJ energy
                    const double sigma2 = sigma * sigma;
                    const double sigma6 = sigma2 * sigma2 * sigma2;
                    const double sigma12 = sigma6 * sigma6;
                    const double r6 = r2_dist * r2_dist * r2_dist;
                    const double r12 = r6 * r6;
                    
                    const double lj_energy = 4.0 * epsilon * (sigma12/r12 - sigma6/r6);
                    
                    platform::log(LogLevel::DEBUG, "LJ pair i=", i, " j=", j, " r=", std::sqrt(r2_dist), " energy=", lj_energy);
                    printf("      Calculated LJ energy = %f kJ/mol (r = %f nm)\n", lj_energy, std::sqrt(r2_dist));
                    
                    // Store in residues
                    if (r1 == r2) {
                        // Intra-residue: all energy to this residue
                        residues[r1].energy_vdw += lj_energy;
                        printf("      Added to residue %d (intra): total vdw = %f\n", 
                               r1, residues[r1].energy_vdw);
                    } else {
                        // Inter-residue: split energy
                        residues[r1].energy_vdw += lj_energy * 0.5;
                        residues[r2].energy_vdw += lj_energy * 0.5;
                        printf("      Split between residues %d and %d (inter)\n", r1, r2);
                    }
                }
            }
        }
    }
}

// Corrected implementation of movement energy using PGP Complete
void computeMovementEnergyPGPCompleteCorrect(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    fprintf(stderr, "\n=== computeMovementEnergyPGPCompleteCorrect called ===\n");
    
    platform::log(LogLevel::DEBUG, "Computing movement energy using corrected PGP Complete method");
    platform::log(LogLevel::DEBUG, "Number of movement residue blocks: ", state.movementResidues.size());
    
    // If no movement residues specified, treat all non-fixed residues as movement
    if (state.movementResidues.empty()) {
        platform::log(LogLevel::DEBUG, "No movement residues specified, treating all non-fixed as movement");
    }
    
    // 1. Calculate grid potential interpolation
    // This gives the interaction of movement atoms with the fixed potential
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate real space electrostatics
    // For movement residues interacting with all residues
    // Skip fixed-fixed pairs (already in grid)
    calculateRealSpacePGPComplete(state, true);
    
    // 3. Calculate self energy for movement atoms only
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, true);
    
    // 4. Calculate LJ for movement residues
    calculateLJPGPComplete(state, true);
    
    // Apply COULOMB constant to real space and residue energies
    state.ewald_energy.real_space *= COULOMB;
    for (auto& res : state.residues) {
        if (res.active) {
            res.energy_elec *= COULOMB;
        }
    }
    
    // Calculate total VDW from movement residues only
    double vdw_total = 0.0;
    printf("\nDEBUG: Calculating total VDW\n");
    printf("  movementResidues.size() = %zu\n", state.movementResidues.size());
    
    if (state.movementResidues.empty()) {
        printf("  No movement residues specified, summing VDW from all non-fixed\n");
        // If no movement residues specified, sum VDW from all non-fixed residues
        for (int i = 0; i < state.activeResidueCount; i++) {
            printf("    Residue %d: active=%d, fixed=%d, energy_vdw=%f\n",
                   i, state.residues[i].active, state.residues[i].fixed, 
                   state.residues[i].energy_vdw);
            if (state.residues[i].active && !state.residues[i].fixed) {
                vdw_total += state.residues[i].energy_vdw;
                printf("      Added to total, now vdw_total = %f\n", vdw_total);
            }
        }
    } else {
        printf("  Using movement residues list\n");
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = movementInfo.startIndex;
                 i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if (state.residues[i].active) {
                    vdw_total += state.residues[i].energy_vdw;
                }
            }
        }
    }
    
    printf("  Final vdw_total = %f\n", vdw_total);
    
    // Set reciprocal (grid) and calculate total
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    platform::log(LogLevel::INFO, "PGP Complete Movement (corrected v2):");
    platform::log(LogLevel::INFO, "  Grid (reciprocal): ", grid_energy);
    platform::log(LogLevel::INFO, "  Real space: ", state.ewald_energy.real_space);
    platform::log(LogLevel::INFO, "  Self: ", state.ewald_energy.self);
    platform::log(LogLevel::INFO, "  VDW: ", vdw_total);
    platform::log(LogLevel::INFO, "  Total: ", state.ewald_energy.total);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc