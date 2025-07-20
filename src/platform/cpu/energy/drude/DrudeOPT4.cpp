/**
 * OPT4 implementation for Drude SCF - 4th order optimization
 * Specifically designed for Drude oscillator model
 */

#include "DrudeForce.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

// Physical constants
static const double ONE_4PI_EPS0 = 138.935456;  // kJ/mol·nm·e^-2

// Forward declaration from OPT3.cpp
bool minimizeDrudePositionsOPT3(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& scfParams,
    const OPT3Coefficients& opt3Coeffs);

/**
 * OPT4 algorithm with 4th order perturbation
 */
static bool minimizeDrudePositionsOPT4(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& scfParams,
    const OPT4Coefficients& opt4Coeffs)
{
    size_t numDrudes = particles.size();
    if (numDrudes == 0) return true;
    
    // Storage for perturbation orders
    std::vector<DrudeForce::Vec3> r0(numDrudes);
    std::vector<DrudeForce::Vec3> r1(numDrudes);
    std::vector<DrudeForce::Vec3> r2(numDrudes);
    std::vector<DrudeForce::Vec3> r3(numDrudes);
    std::vector<DrudeForce::Vec3> r4(numDrudes);
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
        double factor = particles[i].charge / particles[i].kIsotropic;
        r0[i] = electricField[i] * factor;
    }
    
    // First-order
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
    
    // Fourth-order
    std::vector<DrudeForce::Vec3> drudePos3(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos3[i] = parentPos[i] + r3[i];
    }
    calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos3, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        r4[i] = electricField[i] * factor;
    }
    
    // Apply OPT4 formula
    for (size_t i = 0; i < numDrudes; ++i) {
        // Linear combination of perturbation orders
        DrudeForce::Vec3 displacement = 
            r0[i] * opt4Coeffs.c0 + 
            r1[i] * opt4Coeffs.c1 + 
            r2[i] * opt4Coeffs.c2 + 
            r3[i] * opt4Coeffs.c3 +
            r4[i] * opt4Coeffs.c4;
        
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
    
    return true;
}

/**
 * Adaptive OPT algorithm that chooses order based on convergence
 */
static bool minimizeDrudePositionsAdaptiveOPT(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& scfParams,
    const OPT3Coefficients& opt3Coeffs,
    const OPT4Coefficients& opt4Coeffs)
{
    size_t numDrudes = particles.size();
    if (numDrudes == 0) return true;
    
    // Calculate first few orders to determine convergence
    std::vector<DrudeForce::Vec3> r0(numDrudes);
    std::vector<DrudeForce::Vec3> r1(numDrudes);
    std::vector<DrudeForce::Vec3> electricField(numDrudes);
    
    // Get parent positions
    std::vector<DrudeForce::Vec3> parentPos(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        int parentIdx = particles[i].parentIndex;
        parentPos[i] = DrudeForce::Vec3(state.atoms[parentIdx].x,
                                        state.atoms[parentIdx].y,
                                        state.atoms[parentIdx].z);
    }
    
    // Calculate r0 and r1
    calculateElectricFieldAtDrudes(state, particles, screenedPairs, parentPos, electricField, true);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        r0[i] = electricField[i] * factor;
    }
    
    std::vector<DrudeForce::Vec3> drudePos0(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos0[i] = parentPos[i] + r0[i];
    }
    calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos0, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        r1[i] = electricField[i] * factor;
    }
    
    // Analyze convergence
    double avg_convergence_ratio = 0.0;
    int n_valid = 0;
    for (size_t i = 0; i < numDrudes; ++i) {
        double r0_norm = r0[i].norm();
        double r1_norm = r1[i].norm();
        if (r0_norm > 1e-10) {
            avg_convergence_ratio += r1_norm / r0_norm;
            n_valid++;
        }
    }
    
    if (n_valid > 0) {
        avg_convergence_ratio /= n_valid;
    }
    
    // Choose algorithm based on convergence
    if (avg_convergence_ratio < 0.1) {
        // Very fast convergence - use OPT2
        for (size_t i = 0; i < numDrudes; ++i) {
            DrudeForce::Vec3 displacement = r0[i] * 0.2 + r1[i] * 0.8;
            
            double r = displacement.norm();
            if (r > scfParams.maxDrudeDistance) {
                displacement = displacement * (scfParams.maxDrudeDistance / r);
            }
            
            int drudeIdx = particles[i].drudeIndex;
            state.atoms[drudeIdx].x = parentPos[i].x + displacement.x;
            state.atoms[drudeIdx].y = parentPos[i].y + displacement.y;
            state.atoms[drudeIdx].z = parentPos[i].z + displacement.z;
        }
        return true;
    } else if (avg_convergence_ratio < 0.3) {
        // Medium convergence - use OPT3
        return minimizeDrudePositionsOPT3(state, particles, screenedPairs, scfParams, opt3Coeffs);
    } else {
        // Slow convergence - use OPT4
        return minimizeDrudePositionsOPT4(state, particles, screenedPairs, scfParams, opt4Coeffs);
    }
}

/**
 * Hybrid OPT-SCF algorithm
 * For now, just use OPT3 followed by checking displacement magnitude
 */
static bool minimizeDrudePositionsHybrid(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& scfParams,
    const OPT3Coefficients& opt3Coeffs)
{
    // Save initial positions
    std::vector<DrudeForce::Vec3> initialPos(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
        int idx = particles[i].drudeIndex;
        initialPos[i] = DrudeForce::Vec3(state.atoms[idx].x,
                                         state.atoms[idx].y,
                                         state.atoms[idx].z);
    }
    
    // First pass: OPT3 prediction
    minimizeDrudePositionsOPT3(state, particles, screenedPairs, scfParams, opt3Coeffs);
    
    // Check displacement magnitudes
    double maxDisplacement = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        int drudeIdx = particles[i].drudeIndex;
        int parentIdx = particles[i].parentIndex;
        
        DrudeForce::Vec3 drudePos(state.atoms[drudeIdx].x,
                                  state.atoms[drudeIdx].y,
                                  state.atoms[drudeIdx].z);
        DrudeForce::Vec3 parentPos(state.atoms[parentIdx].x,
                                   state.atoms[parentIdx].y,
                                   state.atoms[parentIdx].z);
        
        double displacement = (drudePos - parentPos).norm();
        maxDisplacement = std::max(maxDisplacement, displacement);
    }
    
    // If any displacement is too large, fall back to SCF
    if (maxDisplacement > scfParams.maxDrudeDistance * 0.8) {
        // Restore initial positions and use standard SCF
        for (size_t i = 0; i < particles.size(); ++i) {
            int idx = particles[i].drudeIndex;
            state.atoms[idx].x = initialPos[i].x;
            state.atoms[idx].y = initialPos[i].y;
            state.atoms[idx].z = initialPos[i].z;
        }
        
        // Note: We can't call SCF from here directly
        // Just return false to indicate hybrid failed
        return false;
    }
    
    return true;
}

// Member function implementations
bool DrudeForce::minimizeDrudePositionsWithOPT4(model::MCState& state) {
    return minimizeDrudePositionsOPT4(state, particles, screenedPairs, scfParams, opt4Coeffs);
}

bool DrudeForce::minimizeDrudePositionsAdaptive(model::MCState& state) {
    return minimizeDrudePositionsAdaptiveOPT(state, particles, screenedPairs, scfParams, opt3Coeffs, opt4Coeffs);
}

bool DrudeForce::minimizeDrudePositionsHybrid(model::MCState& state) {
    return ::pygcmc::platform::cpu::minimizeDrudePositionsHybrid(state, particles, screenedPairs, scfParams, opt3Coeffs);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc