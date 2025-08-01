#include "PGPInterpolation.hpp"
#include "PGPGlobal.hpp"
#include "PGPPrecompute.hpp"
#include "platform/platform.hpp"
#include "../common/EnergyUtils.hpp"      // applyPBC
#include "../common/EnergyConstants.hpp"  // COULOMB
#include <cmath>
#include <iostream>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate moving molecule energy through interpolation
 * 
 * This function is another core function of the PGP-PME algorithm, used to quickly evaluate the energy of moving molecules in the precomputed potential field.
 * By using B-spline interpolation from the precomputed grid potential, it avoids direct calculation of intermolecular interactions.
 */
void interpolateMoleculeEnergyImpl(model::MCState& state, double& energy) {
    // Check if parameters are initialized
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Only output logs in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Calculate moving molecule energy through interpolation");
        platform::log(LogLevel::DEBUG, "Number of movement residue groups: ", state.movementResidues.size());
    }
    
    // Reset energy accumulator
    energy = 0.0;
    double raw_energy = 0.0;
    
    // Check if precomputed grid is empty - only perform full check in debug mode
    if (platform::is_debug_mode()) {
        bool gridEmpty = true;
        for (const auto& val : getPGPParams().potentialGrid) {
            if (std::abs(val.real()) > 1e-10 || std::abs(val.imag()) > 1e-10) {
                gridEmpty = false;
                break;
            }
        }
        
        if (gridEmpty) {
            platform::log(LogLevel::WARNING, "PGP grid is empty or not correctly initialized!");
        }
    }
    
    // Process atoms based on movement residue settings
    std::vector<int> residuesToProcess;
    
    if (state.movementResidues.empty()) {
        platform::log(LogLevel::WARNING, "No movement residue information set! Processing all non-fixed residues.");
        
        // Find all non-fixed active residues
        for (size_t i = 0; i < state.residues.size(); i++) {
            if (!state.residues[i].fixed && state.residues[i].active) {
                residuesToProcess.push_back(i);
            }
        }
    } else {
        // Process specified movement residues
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = 0; i < movementInfo.activeCount; i++) {
                int res_idx = movementInfo.startIndex + i;
                if (res_idx < state.activeResidueCount && !state.residues[res_idx].fixed) {
                    residuesToProcess.push_back(res_idx);
                }
            }
        }
    }
    
    // Process atoms in selected residues
    for (int res_idx : residuesToProcess) {
        const auto& residue = state.residues[res_idx];
        
        // Process atoms in this residue
        for (int j = 0; j < residue.atomCount; j++) {
            int atom_index = residue.atomStart + j;
            const auto& atom = state.atoms[atom_index];
            
            // Only process charged atoms
            if (std::abs(atom.charge) < 1e-6) continue;
            
            // Calculate grid position and B-spline interpolation weights
            double pos[3] = {atom.x, atom.y, atom.z};
            
            // Calculate fractional coordinates
            double fractional[3];
            for (int d = 0; d < 3; d++) {
                fractional[d] = pos[d] / getPGPParams().box[d];
                fractional[d] -= floor(fractional[d]);  // Ensure in [0,1) range
                fractional[d] *= getPGPParams().potential_grid_size[d]; // Scale to grid
            }
            
            // Calculate grid index and fractional part
            int gridIndices[3];
            double gridFractions[3];
            for (int d = 0; d < 3; d++) {
                gridFractions[d] = fractional[d] - floor(fractional[d]);
                gridIndices[d] = static_cast<int>(floor(fractional[d]));
                // Ensure grid index within correct range
                if (gridIndices[d] < 0) 
                    gridIndices[d] += getPGPParams().potential_grid_size[d];
            }
            
            // Calculate B-spline coefficients
            int nx = getPGPParams().potential_grid_size[0];
            int ny = getPGPParams().potential_grid_size[1];
            int nz = getPGPParams().potential_grid_size[2];
            int order = getPGPParams().splineOrder;
            
            std::vector<double> thetaX(order);
            std::vector<double> thetaY(order);
            std::vector<double> thetaZ(order);
            
            // Calculate B-spline coefficients for each dimension
            std::vector<double> coefficients(order);
            
            // X dimension B-spline
            computeBSplineCoefficients(gridFractions[0], order, coefficients);
            for (int k = 0; k < order; k++) {
                thetaX[k] = coefficients[k];
            }
            
            // Y dimension B-spline
            computeBSplineCoefficients(gridFractions[1], order, coefficients);
            for (int k = 0; k < order; k++) {
                thetaY[k] = coefficients[k];
            }
            
            // Z dimension B-spline
            computeBSplineCoefficients(gridFractions[2], order, coefficients);
            for (int k = 0; k < order; k++) {
                thetaZ[k] = coefficients[k];
            }
            
            // Interpolate potential
            double potential = 0.0;
            
            // Loop through all B-spline support points
            for (int ix = 0; ix < order; ix++) {
                int xindex = (gridIndices[0] + ix) % nx;
                
                for (int iy = 0; iy < order; iy++) {
                    int yindex = (gridIndices[1] + iy) % ny;
                    
                    for (int iz = 0; iz < order; iz++) {
                        int zindex = (gridIndices[2] + iz) % nz;
                        
                        // Calculate three-dimensional grid index
                        int index = xindex * ny * nz + yindex * nz + zindex;
                        
                        // Use B-spline weights to accumulate potential
                        double grid_value = getPGPParams().potentialGrid[index].real();
                        double weight = thetaX[ix] * thetaY[iy] * thetaZ[iz];
                        potential += grid_value * weight;
                    }
                }
            }
            
            // Accumulate energy (potential * charge)
            double atom_energy = potential * atom.charge;
            raw_energy += atom_energy;
            
            if (platform::is_debug_mode()) {
                platform::log(LogLevel::DEBUG, "Atom potential: ", potential, ", Atom energy contribution: ", atom_energy);
            }
        }
    }
    
    // The potential already includes the correct normalization factor from precomputation
    energy = raw_energy;
    
    // Only output energy calculation details in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Final calculated PGP energy: ", energy, " kJ/mol");
    }
}

/**
 * @brief Calculate moving molecule energy through interpolation and return calculation result
 */
double calculateMoleculeEnergyImpl(model::MCState& state) {
    // Reset residue energies to avoid carry-over between moves
    for (int r = 0; r < state.activeResidueCount; ++r) {
        if (state.residues[r].active) {
            state.residues[r].energy_elec = 0.0f;
            state.residues[r].energy_vdw = 0.0f;
        }
    }
    
    double energy = 0.0;
    interpolateMoleculeEnergyImpl(state, energy);
    
    /* ---------------- add mov-mov reciprocal correction ---------------- */
    // Collect atom indices of movement residues with residue info
    struct AtomInfo {
        int atomIdx;
        int resIdx;
    };
    std::vector<AtomInfo> movAtoms;
    
    auto isMovementResidue = [&](int resIdx) {
        if (state.movementResidues.empty())          // treat all non-fixed as movement
            return !state.residues[resIdx].fixed;
        for (const auto& m : state.movementResidues)
            if (resIdx >= m.startIndex && resIdx < m.startIndex + m.activeCount) return true;
        return false;
    };

    for (int r = 0; r < state.activeResidueCount; ++r) {
        if (!isMovementResidue(r) || !state.residues[r].active) continue;
        for (int i = 0; i < state.residues[r].atomCount; ++i) {
            AtomInfo info;
            info.atomIdx = state.residues[r].atomStart + i;
            info.resIdx = r;
            movAtoms.push_back(info);
        }
    }

    const double alpha = getPGPParams().alpha;
    const float* box   = state.info.box;
    const auto& atoms  = state.atoms;

    double movMovRecip = 0.0;
    for (size_t a = 0; a + 1 < movAtoms.size(); ++a) {
        int i = movAtoms[a].atomIdx;
        int resI = movAtoms[a].resIdx;
        const double qi = atoms[i].charge;
        if (qi == 0.0) continue;

        for (size_t b = a + 1; b < movAtoms.size(); ++b) {
            int j = movAtoms[b].atomIdx;
            int resJ = movAtoms[b].resIdx;
            
            // Skip intramolecular interactions
            if (resI == resJ) continue;
            
            const double qj = atoms[j].charge;
            if (qj == 0.0) continue;

            double dx = atoms[i].x - atoms[j].x;
            double dy = atoms[i].y - atoms[j].y;
            double dz = atoms[i].z - atoms[j].z;

            float dx_f = static_cast<float>(dx);
            float dy_f = static_cast<float>(dy);
            float dz_f = static_cast<float>(dz);
            applyPBC(dx_f, dy_f, dz_f, box);
            dx = dx_f; dy = dy_f; dz = dz_f;

            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;

            const double r      = std::sqrt(r2);
            const double energy_pair = qi * qj * std::erf(alpha * r) / r;
            movMovRecip += energy_pair;
        }
    }
    energy += COULOMB * movMovRecip;
    /* ------------------------------------------------------------------- */
    
    // If no fixed residues, potential grid may need to be recomputed
    if (std::abs(energy) < 1e-10) {
        int fixed_count = 0;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].fixed && state.residues[i].active) {
                fixed_count++;
            }
        }
        
        if (fixed_count == 0) {
            platform::log(LogLevel::WARNING, "No fixed residues found, energy near zero!");
            // Call the public interface instead of the impl to avoid lock issues
            precomputeGridPotential(state, false);
            interpolateMoleculeEnergyImpl(state, energy);
        }
    }
    
    return energy;
}

double computeMoleculeEnergyGlobalImpl(model::MCState& state, const std::vector<int>& movementResidues, 
                                       const std::vector<int>& nearbyResidues, int threadIndex) {
    // Suppress unused parameter warnings
    (void)nearbyResidues;
    (void)threadIndex;
    
    // If parameters are not initialized, return 0
    if (!getPGPParams().initialized) {
        platform::log(LogLevel::WARNING, "PGP parameters not initialized, returning 0 energy");
        return 0.0;
    }
    
    double totalEnergy = 0.0;
    
    // Process different residue energy calculations
    if (movementResidues.empty()) {
        totalEnergy = calculateMoleculeEnergyImpl(state);
    } else {
        // If specified movement residues, we need to modify state's movementResidues
        auto originalMovementResidues = state.movementResidues;
        
        // Clear and set new movementResidues
        state.movementResidues.clear();
        model::MCMovementResidueInfo info;
        info.startIndex = movementResidues[0];
        info.activeCount = movementResidues.size();
        state.movementResidues.push_back(info);
        
        // Calculate energy
        totalEnergy = calculateMoleculeEnergyImpl(state);
        
        // Restore original movementResidues
        state.movementResidues = originalMovementResidues;
    }
    
    return totalEnergy;
}

// === Public Interface Wrappers ===

void interpolateMoleculeEnergy(model::MCState& state, double& energy) {
    interpolateMoleculeEnergyImpl(state, energy);
}

double calculateMoleculeEnergy(model::MCState& state) {
    return calculateMoleculeEnergyImpl(state);
}

double computeMoleculeEnergyGlobal(model::MCState& state, const std::vector<int>& movementResidues, 
                                   const std::vector<int>& nearbyResidues, int threadIndex) {
    return computeMoleculeEnergyGlobalImpl(state, movementResidues, nearbyResidues, threadIndex);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 