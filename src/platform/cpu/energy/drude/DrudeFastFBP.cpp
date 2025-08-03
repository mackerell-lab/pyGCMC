/**
 * @file DrudeFastFBP.cpp
 * @brief Fast Force Balance Predictor implementation
 */

#include "DrudeFastFBP.hpp"
#include "../common/EnergyConstants.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

bool DrudeFastFBP::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    // Mark screenedPairs as unused for now (can be implemented later)
    (void)screenedPairs;
    
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};
    
    // Initialize statistics
    stats_ = {0, 0.0, 0.0, false};
    
    // Configure for system if adaptive mode enabled
    if (adaptiveMode_) {
        configureForSystem(state, particles.size());
    }
    
    // Store positions for convergence check
    std::vector<Vec3> previousPositions(particles.size());
    std::vector<double> errorHistory;
    
    // Determine number of iterations based on mode
    int maxIter = (iterMode_ == IterationMode::Fixed) ? fbpIterations_ : maxIterations_;
    
    // Main FBP iterations
    for (int iter = 0; iter < maxIter; ++iter) {
        
        // Save previous positions
        saveDrudePositions(state, particles, previousPositions);
        
        // Update each Drude position
        for (size_t i = 0; i < particles.size(); ++i) {
            const auto& particle = particles[i];
            const auto& parent = state.atoms[particle.parentIndex];
            auto& drude = state.atoms[particle.drudeIndex];
            
            Vec3 parentPos = {parent.x, parent.y, parent.z};
            Vec3 drudePos = {drude.x, drude.y, drude.z};
            
            // Step 1: Compute field from fixed charges at Drude position
            Vec3 fieldFixed = computeFixedField(state, drudePos, 
                                              particle.drudeIndex, 
                                              particle.parentIndex);
            
            // Step 2: Compute field from other Drude particles (if enabled)
            Vec3 fieldDrude = {0.0, 0.0, 0.0};
            if (includeDrudeDrude_ && iter > 0) {
                fieldDrude = computeDrudeField(state, particles, drudePos, i);
            }
            
            // Step 3: Total field
            Vec3 fieldTotal = {
                fieldFixed[0] + fieldDrude[0],
                fieldFixed[1] + fieldDrude[1],
                fieldFixed[2] + fieldDrude[2]
            };
            
            // Step 4: Force balance: F_spring + F_electric = 0
            // For equilibrium: k*(r_drude - r_parent) + q*E = 0
            // This gives: displacement = (q/k)*E
            // Note: For negative charge, displacement is opposite to field
            double dispFactor = particle.charge / particle.kSpring;
            
            Vec3 displacement = {
                dispFactor * fieldTotal[0],
                dispFactor * fieldTotal[1],
                dispFactor * fieldTotal[2]
            };
            
            // Apply damping for stability
            Vec3 currentDisp = {
                drudePos[0] - parentPos[0],
                drudePos[1] - parentPos[1],
                drudePos[2] - parentPos[2]
            };
            
            // Use adaptive damping if enabled
            double effectiveDamping = (iter == 0) ? 1.0 : 
                (adaptiveMode_ ? adaptiveParams_.dampingFactor : dampingFactor_);
            
            Vec3 newDisp = {
                effectiveDamping * displacement[0] + (1 - effectiveDamping) * currentDisp[0],
                effectiveDamping * displacement[1] + (1 - effectiveDamping) * currentDisp[1],
                effectiveDamping * displacement[2] + (1 - effectiveDamping) * currentDisp[2]
            };
            
            Vec3 newDrudePos = {
                parentPos[0] + newDisp[0],
                parentPos[1] + newDisp[1],
                parentPos[2] + newDisp[2]
            };
            
            // Step 5: Apply constraints if needed
            if (params.enableHardWall && params.maxDrudeDistance > 0) {
                applyConstraints(newDrudePos, parentPos, params.maxDrudeDistance);
            }
            
            // Update position
            drude.x = newDrudePos[0];
            drude.y = newDrudePos[1];
            drude.z = newDrudePos[2];
        }
        
        // Check convergence for dynamic mode
        if (iterMode_ != IterationMode::Fixed && iter >= minIterations_) {
            double forceTolerance = 0.1; // kJ/mol/nm
            bool converged = checkConvergence(state, particles, previousPositions, 
                                            convTolerance_, forceTolerance);
            
            if (converged) {
                stats_.actualIterations = iter + 1;
                stats_.converged = true;
                stats_.finalError = convTolerance_;
                return true;
            }
            
            // Calculate convergence rate for adaptive mode
            if (iterMode_ == IterationMode::Adaptive && iter > minIterations_) {
                double maxDisp = 0.0;
                for (size_t i = 0; i < particles.size(); ++i) {
                    const auto& drude = state.atoms[particles[i].drudeIndex];
                    Vec3 newPos = {drude.x, drude.y, drude.z};
                    
                    double dx = newPos[0] - previousPositions[i][0];
                    double dy = newPos[1] - previousPositions[i][1];
                    double dz = newPos[2] - previousPositions[i][2];
                    
                    // Apply PBC to displacement
                    if (dx > halfBox[0]) dx -= box[0];
                    if (dx < -halfBox[0]) dx += box[0];
                    if (dy > halfBox[1]) dy -= box[1];
                    if (dy < -halfBox[1]) dy += box[1];
                    if (dz > halfBox[2]) dz -= box[2];
                    if (dz < -halfBox[2]) dz += box[2];
                    
                    double disp2 = dx*dx + dy*dy + dz*dz;
                    maxDisp = std::max(maxDisp, std::sqrt(disp2));
                }
                
                errorHistory.push_back(maxDisp);
                
                // Adjust adaptive parameters
                if (errorHistory.size() >= 3) {
                    double convergenceRate = calculateConvergenceRate(errorHistory);
                    adjustAdaptiveParameters(convergenceRate);
                    stats_.convergenceRate = convergenceRate;
                }
            }
        }
        
        // Update statistics
        stats_.actualIterations = iter + 1;
    }
    
    // Final statistics
    stats_.converged = (iterMode_ == IterationMode::Fixed); // Fixed mode always "converges"
    return true;
}

Vec3 DrudeFastFBP::computeFixedField(
    const model::MCState& state,
    const Vec3& drudePos,
    int drudeIdx,
    int parentIdx
) {
    Vec3 field = {0.0, 0.0, 0.0};
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};
    
    // Sum over all atoms except self
    for (int i = 0; i < state.activeAtomCount; ++i) {
        // Skip self only (parent contributes to field)
        if (i == drudeIdx) continue;
        
        const auto& atom = state.atoms[i];
        
        // Skip if no charge
        if (std::abs(atom.charge) < 1e-6) continue;
        
        // Skip parent-Drude interaction (handled by spring)
        if (i == parentIdx) continue;
        
        // Skip intramolecular interactions
        // Drude particles should NOT interact with other atoms in the same molecule
        if (inSameMolecule(drudeIdx, i, state)) continue;
        
        // Compute distance with PBC
        double dx = drudePos[0] - atom.x;
        double dy = drudePos[1] - atom.y;
        double dz = drudePos[2] - atom.z;
        
        // Apply minimum image convention
        if (dx > halfBox[0]) dx -= box[0];
        if (dx < -halfBox[0]) dx += box[0];
        if (dy > halfBox[1]) dy -= box[1];
        if (dy < -halfBox[1]) dy += box[1];
        if (dz > halfBox[2]) dz -= box[2];
        if (dz < -halfBox[2]) dz += box[2];
        
        double r2 = dx*dx + dy*dy + dz*dz;
        
        // Apply cutoff
        if (r2 > cutoff2 || r2 < 1e-10) continue;
        
        // E = k*q/r^2 * r_hat
        double r = std::sqrt(r2);
        double fieldMag = DrudeConstants::ONE_4PI_EPS0 * atom.charge / (r2 * r);
        
        field[0] += fieldMag * dx;
        field[1] += fieldMag * dy;
        field[2] += fieldMag * dz;
    }
    
    return field;
}

Vec3 DrudeFastFBP::computeDrudeField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const Vec3& drudePos,
    int currentDrudeIdx
) {
    Vec3 field = {0.0, 0.0, 0.0};
    const double cutoff2 = drudeCutoff_ * drudeCutoff_;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};
    
    // Sum over other Drude particles
    for (size_t i = 0; i < particles.size(); ++i) {
        if (static_cast<int>(i) == currentDrudeIdx) continue;
        
        const auto& otherDrude = state.atoms[particles[i].drudeIndex];
        
        // Compute distance with PBC
        double dx = drudePos[0] - otherDrude.x;
        double dy = drudePos[1] - otherDrude.y;
        double dz = drudePos[2] - otherDrude.z;
        
        // Apply minimum image convention
        if (dx > halfBox[0]) dx -= box[0];
        if (dx < -halfBox[0]) dx += box[0];
        if (dy > halfBox[1]) dy -= box[1];
        if (dy < -halfBox[1]) dy += box[1];
        if (dz > halfBox[2]) dz -= box[2];
        if (dz < -halfBox[2]) dz += box[2];
        
        double r2 = dx*dx + dy*dy + dz*dz;
        
        // Apply cutoff for Drude-Drude
        if (r2 > cutoff2 || r2 < 1e-10) continue;
        
        // E = k*q/r^2 * r_hat
        double r = std::sqrt(r2);
        double fieldMag = DrudeConstants::ONE_4PI_EPS0 * particles[i].charge / (r2 * r);
        
        field[0] += fieldMag * dx;
        field[1] += fieldMag * dy;
        field[2] += fieldMag * dz;
    }
    
    return field;
}

void DrudeFastFBP::applyConstraints(
    Vec3& drudePos,
    const Vec3& parentPos,
    double maxDistance
) {
    // Calculate displacement
    double dx = drudePos[0] - parentPos[0];
    double dy = drudePos[1] - parentPos[1];
    double dz = drudePos[2] - parentPos[2];
    
    double r2 = dx*dx + dy*dy + dz*dz;
    double maxDist2 = maxDistance * maxDistance;
    
    // If within constraint, no action needed
    if (r2 <= maxDist2) {
        return;
    }
    
    // Scale back to maximum distance
    double r = std::sqrt(r2);
    double scale = maxDistance / r;
    
    drudePos[0] = parentPos[0] + scale * dx;
    drudePos[1] = parentPos[1] + scale * dy;
    drudePos[2] = parentPos[2] + scale * dz;
}

bool DrudeFastFBP::inSameMolecule(
    int atom1,
    int atom2,
    const model::MCState& state
) {
    // Find which residue each atom belongs to
    int res1 = -1, res2 = -1;
    
    for (int i = 0; i < state.activeResidueCount; ++i) {
        const auto& res = state.residues[i];
        if (atom1 >= res.atomStart && atom1 < res.atomStart + res.atomCount) {
            res1 = i;
        }
        if (atom2 >= res.atomStart && atom2 < res.atomStart + res.atomCount) {
            res2 = i;
        }
    }
    
    return (res1 >= 0 && res1 == res2);
}

bool DrudeFastFBP::checkConvergence(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& previousPositions,
    double dispTolerance,
    double forceTolerance
) {
    double maxDisplacement = 0.0;
    double maxForceImbalance = 0.0;
    
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& drude = state.atoms[particle.drudeIndex];
        
        // Check displacement change
        double dx = drude.x - previousPositions[i][0];
        double dy = drude.y - previousPositions[i][1];
        double dz = drude.z - previousPositions[i][2];
        
        // Apply PBC to displacement
        if (dx > halfBox[0]) dx -= box[0];
        if (dx < -halfBox[0]) dx += box[0];
        if (dy > halfBox[1]) dy -= box[1];
        if (dy < -halfBox[1]) dy += box[1];
        if (dz > halfBox[2]) dz -= box[2];
        if (dz < -halfBox[2]) dz += box[2];
        
        double displacement = std::sqrt(dx*dx + dy*dy + dz*dz);
        maxDisplacement = std::max(maxDisplacement, displacement);
        
        // Check force balance
        Vec3 netForce = calculateNetForce(state, particle, i, particles);
        double forceMag = std::sqrt(
            netForce[0]*netForce[0] + 
            netForce[1]*netForce[1] + 
            netForce[2]*netForce[2]
        );
        maxForceImbalance = std::max(maxForceImbalance, forceMag);
    }
    
    // Double convergence criteria
    bool dispConverged = maxDisplacement < dispTolerance;
    bool forceConverged = maxForceImbalance < forceTolerance;
    
    return dispConverged && forceConverged;
}

void DrudeFastFBP::saveDrudePositions(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& positions
) {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& drude = state.atoms[particles[i].drudeIndex];
        positions[i] = {drude.x, drude.y, drude.z};
    }
}

Vec3 DrudeFastFBP::calculateNetForce(
    const model::MCState& state,
    const DrudeParticle& particle,
    size_t particleIndex,
    const std::vector<DrudeParticle>& allParticles
) {
    const auto& drude = state.atoms[particle.drudeIndex];
    const auto& parent = state.atoms[particle.parentIndex];
    
    Vec3 drudePos = {drude.x, drude.y, drude.z};
    Vec3 parentPos = {parent.x, parent.y, parent.z};
    
    // Spring force
    Vec3 springForce = {
        particle.kSpring * (parentPos[0] - drudePos[0]),
        particle.kSpring * (parentPos[1] - drudePos[1]),
        particle.kSpring * (parentPos[2] - drudePos[2])
    };
    
    // Electric field at Drude position
    Vec3 fieldFixed = computeFixedField(state, drudePos, 
                                      particle.drudeIndex, 
                                      particle.parentIndex);
    
    Vec3 fieldDrude = computeDrudeField(state, allParticles, drudePos, particleIndex);
    
    // Electric force
    Vec3 electricForce = {
        particle.charge * (fieldFixed[0] + fieldDrude[0]),
        particle.charge * (fieldFixed[1] + fieldDrude[1]),
        particle.charge * (fieldFixed[2] + fieldDrude[2])
    };
    
    // Net force
    return {
        springForce[0] + electricForce[0],
        springForce[1] + electricForce[1],
        springForce[2] + electricForce[2]
    };
}

void DrudeFastFBP::configureForSystem(
    const model::MCState& state,
    size_t nParticles
) {
    // Calculate system density
    double volume = state.info.box[0] * state.info.box[1] * state.info.box[2];
    double density = nParticles / volume;
    
    // Adjust parameters based on system size and density
    if (nParticles < 10) {
        // Small system: use fast settings
        fbpIterations_ = 3;
        drudeCutoff_ = 0.8;
        dampingFactor_ = 0.7;
        adaptiveParams_.dampingFactor = 0.7;
        adaptiveParams_.cutoff = 0.8;
    } else if (nParticles < 50) {
        // Medium system: balanced settings
        fbpIterations_ = 5;
        drudeCutoff_ = 1.0;
        dampingFactor_ = 0.75;
        adaptiveParams_.dampingFactor = 0.75;
        adaptiveParams_.cutoff = 1.0;
    } else {
        // Large system: use system cutoff
        fbpIterations_ = 7;
        drudeCutoff_ = std::min(static_cast<double>(state.info.cutoff), 1.4);
        dampingFactor_ = 0.8;
        adaptiveParams_.dampingFactor = 0.8;
        adaptiveParams_.cutoff = drudeCutoff_;
    }
    
    // Adjust for high density
    if (density > 30.0) { // High density threshold
        fbpIterations_ += 2;
        dampingFactor_ = std::min(dampingFactor_ + 0.05, 0.9);
    }
}

void DrudeFastFBP::adjustAdaptiveParameters(
    double convergenceRate
) {
    adaptiveParams_.convergenceRate = convergenceRate;
    
    if (convergenceRate < 0.1) {
        // Convergence too slow, reduce damping
        adaptiveParams_.dampingFactor *= 0.9;
        adaptiveParams_.stagnationCount++;
        
        // If stagnating, increase cutoff
        if (adaptiveParams_.stagnationCount > 3) {
            adaptiveParams_.cutoff = std::min(adaptiveParams_.cutoff + 0.1, 1.4);
            adaptiveParams_.stagnationCount = 0;
        }
    } else if (convergenceRate > 0.5) {
        // Convergence possibly oscillating, increase damping
        adaptiveParams_.dampingFactor *= 1.1;
        adaptiveParams_.dampingFactor = std::min(adaptiveParams_.dampingFactor, 0.95);
        adaptiveParams_.stagnationCount = 0;
    } else {
        // Good convergence rate
        adaptiveParams_.stagnationCount = 0;
    }
}

double DrudeFastFBP::calculateConvergenceRate(
    const std::vector<double>& errorHistory
) {
    int n = errorHistory.size();
    if (n < 3) return 1.0;
    
    // Calculate improvement rate from recent iterations
    double improvement1 = (errorHistory[n-2] - errorHistory[n-1]) / errorHistory[n-2];
    double improvement2 = (errorHistory[n-3] - errorHistory[n-2]) / errorHistory[n-3];
    
    return (improvement1 + improvement2) / 2.0;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc