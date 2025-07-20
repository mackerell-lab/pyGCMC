/**
 * OPT3 implementation for Drude SCF
 * Based on "An empirical extrapolation scheme for efficient treatment of induced dipoles"
 * Adapted for Drude oscillator model
 */

#include "DrudeForce.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

// Physical constants
static const double ONE_4PI_EPS0 = 138.935456;  // kJ/mol·nm·e^-2

/**
 * Calculate electric field at Drude positions including Thole screening
 * This is the core operation needed for each perturbation order
 */
void calculateElectricFieldAtDrudes(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const std::vector<DrudeForce::Vec3>& drudePositions,
    std::vector<DrudeForce::Vec3>& electricField,
    bool includeStaticOnly)
{
    // Clear electric field
    std::fill(electricField.begin(), electricField.end(), DrudeForce::Vec3());
    
    // Calculate field from all charges
    for (size_t i = 0; i < particles.size(); ++i) {
        int drudeIdx = particles[i].drudeIndex;
        DrudeForce::Vec3 posDrude = drudePositions[i];
        
        // Field from all atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self-interaction
            if (j == drudeIdx) continue;
            
            // Skip atoms in the same residue (intramolecular exclusion)
            int drudeParentIdx = particles[i].parentIndex;
            bool sameResidue = false;
            for (int resIdx = 0; resIdx < state.activeResidueCount; ++resIdx) {
                const auto& res = state.residues[resIdx];
                bool drudeInRes = (drudeParentIdx >= res.atomStart && 
                                   drudeParentIdx < res.atomStart + res.atomCount);
                bool atomInRes = (j >= res.atomStart && 
                                 j < res.atomStart + res.atomCount);
                if (drudeInRes && atomInRes) {
                    sameResidue = true;
                    break;
                }
            }
            if (sameResidue) continue;
            
            // For static-only, skip other Drude particles
            if (includeStaticOnly) {
                bool isDrude = false;
                for (const auto& p : particles) {
                    if (j == p.drudeIndex) {
                        isDrude = true;
                        break;
                    }
                }
                if (isDrude) continue;
            }
            
            DrudeForce::Vec3 posJ(state.atoms[j].x, state.atoms[j].y, state.atoms[j].z);
            DrudeForce::Vec3 delta = posJ - posDrude;
            
            // Apply PBC if needed
            if (state.info.box[0] > 0) {
                delta.x -= state.info.box[0] * std::round(delta.x / state.info.box[0]);
                delta.y -= state.info.box[1] * std::round(delta.y / state.info.box[1]);
                delta.z -= state.info.box[2] * std::round(delta.z / state.info.box[2]);
            }
            
            double r = delta.norm();
            if (r < 1e-10) continue;
            
            // Check cutoff
            if (state.info.cutoff > 0 && r > state.info.cutoff) continue;
            
            // Check if this is a screened interaction
            bool isScreened = false;
            double thole = 0.0;
            
            // Find if this pair is screened
            for (const auto& pair : screenedPairs) {
                // Check if this is a screened Drude-Drude interaction
                if ((i == static_cast<size_t>(pair.dipole1) && j == particles[pair.dipole2].drudeIndex) ||
                    (i == static_cast<size_t>(pair.dipole2) && j == particles[pair.dipole1].drudeIndex)) {
                    isScreened = true;
                    thole = pair.thole;
                    break;
                }
            }
            
            // Calculate field with or without screening
            if (isScreened && !includeStaticOnly) {
                // Thole screening for Drude-Drude interactions
                double alpha_i = particles[i].polarizability;
                double alpha_j = 0.0;
                // Find alpha_j
                for (size_t k = 0; k < particles.size(); ++k) {
                    if (particles[k].drudeIndex == j) {
                        alpha_j = particles[k].polarizability;
                        break;
                    }
                }
                
                if (alpha_j > 0) {
                    double u = r / std::pow(alpha_i * alpha_j, 1.0/6.0);
                    double screening = 1.0 - (1.0 + thole*u/2.0) * std::exp(-thole*u);
                    double fieldMag = ONE_4PI_EPS0 * state.atoms[j].charge * screening / (r * r);
                    electricField[i] += delta * (fieldMag / r);
                } else {
                    // Regular interaction if not Drude
                    double fieldMag = ONE_4PI_EPS0 * state.atoms[j].charge / (r * r);
                    electricField[i] += delta * (fieldMag / r);
                }
            } else {
                // Regular Coulomb field
                double fieldMag = ONE_4PI_EPS0 * state.atoms[j].charge / (r * r);
                electricField[i] += delta * (fieldMag / r);
            }
        }
    }
}

/**
 * OPT3 algorithm for Drude position optimization
 * Returns true if successful
 */
bool minimizeDrudePositionsOPT3(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& scfParams,
    const OPT3Coefficients& opt3Coeffs)
{
    const OPT3Coefficients& coeff = opt3Coeffs;
    size_t numDrudes = particles.size();
    
    // Storage for different order displacements
    std::vector<DrudeForce::Vec3> r0(numDrudes);  // Zero-order
    std::vector<DrudeForce::Vec3> r1(numDrudes);  // First-order
    std::vector<DrudeForce::Vec3> r2(numDrudes);  // Second-order
    std::vector<DrudeForce::Vec3> r3(numDrudes);  // Third-order
    std::vector<DrudeForce::Vec3> electricField(numDrudes);
    
    // Get parent positions
    std::vector<DrudeForce::Vec3> parentPos(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        int parentIdx = particles[i].parentIndex;
        parentPos[i] = DrudeForce::Vec3(state.atoms[parentIdx].x,
                                        state.atoms[parentIdx].y,
                                        state.atoms[parentIdx].z);
    }
    
    // Zero-order: response to static field only
    calculateElectricFieldAtDrudes(state, particles, screenedPairs, parentPos, electricField, true);
    for (size_t i = 0; i < numDrudes; ++i) {
        // r0 = (q_D/k) * E_static
        double factor = particles[i].charge / particles[i].kIsotropic;
        r0[i] = electricField[i] * factor;
    }
    
    // First-order: include field from zero-order Drudes
    std::vector<DrudeForce::Vec3> drudePos0(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos0[i] = parentPos[i] + r0[i];
    }
    calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos0, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        r1[i] = electricField[i] * factor;
    }
    
    // Second-order
    std::vector<DrudeForce::Vec3> drudePos1(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos1[i] = parentPos[i] + r1[i];
    }
    calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos1, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        r2[i] = electricField[i] * factor;
    }
    
    // Third-order
    std::vector<DrudeForce::Vec3> drudePos2(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos2[i] = parentPos[i] + r2[i];
    }
    calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos2, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        r3[i] = electricField[i] * factor;
    }
    
    // Combine using OPT3 coefficients
    for (size_t i = 0; i < numDrudes; ++i) {
        DrudeForce::Vec3 displacement = r0[i] * coeff.c0 + 
                                       r1[i] * coeff.c1 + 
                                       r2[i] * coeff.c2 + 
                                       r3[i] * coeff.c3;
        
        // Apply hard wall constraint
        double r = displacement.norm();
        if (r > scfParams.maxDrudeDistance) {
            displacement = displacement * (scfParams.maxDrudeDistance / r);
        }
        
        // Update Drude position
        int drudeIdx = particles[i].drudeIndex;
        state.atoms[drudeIdx].x = parentPos[i].x + displacement.x;
        state.atoms[drudeIdx].y = parentPos[i].y + displacement.y;
        state.atoms[drudeIdx].z = parentPos[i].z + displacement.z;
    }
    
    return true;  // OPT3 always "converges"
}

// Add this as an alternative minimization method in DrudeForce
bool DrudeForce::minimizeDrudePositionsWithOPT3(model::MCState& state) {
    return minimizeDrudePositionsOPT3(state, particles, screenedPairs, scfParams, opt3Coeffs);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc