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
    
    // Store initial positions for convergence check
    std::vector<Vec3> oldPositions(particles.size());
    
    // Main FBP iterations
    for (int iter = 0; iter < fbpIterations_; ++iter) {
        
        // Save old positions
        for (size_t i = 0; i < particles.size(); ++i) {
            const auto& drude = state.atoms[particles[i].drudeIndex];
            oldPositions[i] = {drude.x, drude.y, drude.z};
        }
        
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
            
            // Apply damping for stability (use class member)
            Vec3 currentDisp = {
                drudePos[0] - parentPos[0],
                drudePos[1] - parentPos[1],
                drudePos[2] - parentPos[2]
            };
            
            // First iteration: direct placement, later: damped update
            double effectiveDamping = (iter == 0) ? 1.0 : dampingFactor_;
            
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
        
        // Check convergence (only on last iteration for speed)
        if (iter == fbpIterations_ - 1) {
            double maxDisp = 0.0;
            for (size_t i = 0; i < particles.size(); ++i) {
                const auto& drude = state.atoms[particles[i].drudeIndex];
                Vec3 newPos = {drude.x, drude.y, drude.z};
                
                double dx = newPos[0] - oldPositions[i][0];
                double dy = newPos[1] - oldPositions[i][1];
                double dz = newPos[2] - oldPositions[i][2];
                
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
            
            // Very loose convergence for speed (5% accuracy target)
            const double looseTol = 0.001; // 0.001 nm = 0.01 Å
            if (maxDisp < looseTol) {
                return true; // Converged
            }
        }
    }
    
    return true; // Always return true for FastFBP (best effort)
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

} // namespace cpu
} // namespace platform
} // namespace pygcmc