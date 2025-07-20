// DrudeForceBalance.cpp - Force Balance Predictor for Drude particles

#include "DrudeForce.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * Force Balance Predictor (FBP) - Direct solution based on force equilibrium
 * 
 * Key insight: At equilibrium, F_spring + F_electric + F_thole = 0
 * This gives: r_drude = r_parent + (q*E + F_thole) / k
 */
bool DrudeForce::minimizeDrudePositionsWithFBP(model::MCState& state) {
    const int maxIterations = scfParams.maxIterations;  // Use SCF parameters
    const double convergenceTolerance = scfParams.tolerance;  // Use SCF tolerance
    double dampingFactor = scfParams.dampingFactor;  // Use SCF damping
    
    // Store current positions for convergence check
    std::vector<Vec3> oldPositions(particles.size());
    std::vector<Vec3> drudeDisplacements(particles.size());
    
    // Step 1: Calculate fixed electric field (from non-Drude atoms)
    std::vector<Vec3> fixedField(particles.size());
    calculateFixedElectricField(state, fixedField);
    
    // Step 2: Initial guess - ignore Drude-Drude interactions
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        double k = particle.kIsotropic;  // Use actual force constant
        
        // Initial displacement from fixed field only
        // F = q*E, displacement = F/k = q*E/k
        Vec3 displacement = fixedField[i] * (particle.charge / k);
        
        // Update Drude position
        state.atoms[particle.drudeIndex].x = state.atoms[particle.parentIndex].x + displacement.x;
        state.atoms[particle.drudeIndex].y = state.atoms[particle.parentIndex].y + displacement.y;
        state.atoms[particle.drudeIndex].z = state.atoms[particle.parentIndex].z + displacement.z;
        
        drudeDisplacements[i] = displacement;
    }
    
    // Step 3: Iterative refinement with Drude-Drude interactions
    bool converged = false;
    double prevMaxForce = 1e10;
    int stuckCount = 0;
    
    for (int iter = 0; iter < maxIterations; iter++) {
        // Save old positions
        for (size_t i = 0; i < particles.size(); i++) {
            const auto& drude = state.atoms[particles[i].drudeIndex];
            oldPositions[i] = Vec3(drude.x, drude.y, drude.z);
        }
        
        // Calculate total forces on each Drude
        std::vector<Vec3> totalForces(particles.size());
        calculateDrudeForces(state, totalForces, fixedField);
        
        // Update positions based on force balance
        double maxForce = 0.0;
        
        for (size_t i = 0; i < particles.size(); i++) {
            const auto& particle = particles[i];
            double k = particle.kIsotropic;  // Use actual force constant
            
            // Force balance: F_spring + F_external = 0
            // -k*(r_d - r_p) + F_external = 0
            // r_d = r_p + F_external/k
            
            Vec3 newDisplacement = totalForces[i] * (1.0 / k);
            
            // Apply damping for stability
            Vec3 dampedDisplacement = drudeDisplacements[i] * (1.0 - dampingFactor) + 
                                     newDisplacement * dampingFactor;
            
            // Limit displacement magnitude (hard wall constraint)
            double dispMag = dampedDisplacement.norm();
            if (dispMag > scfParams.maxDrudeDistance) {
                dampedDisplacement = dampedDisplacement * (scfParams.maxDrudeDistance / dispMag);
            }
            
            // Update position
            state.atoms[particle.drudeIndex].x = state.atoms[particle.parentIndex].x + dampedDisplacement.x;
            state.atoms[particle.drudeIndex].y = state.atoms[particle.parentIndex].y + dampedDisplacement.y;
            state.atoms[particle.drudeIndex].z = state.atoms[particle.parentIndex].z + dampedDisplacement.z;
            
            drudeDisplacements[i] = dampedDisplacement;
            
            // Check convergence - calculate actual force on Drude
            // F_total = F_external - k*(r_d - r_p)
            Vec3 springForce = dampedDisplacement * k;
            Vec3 residualForce = totalForces[i] - springForce;
            maxForce = std::max(maxForce, residualForce.norm());
        }
        
        // Check convergence
        if (maxForce < convergenceTolerance) {
            converged = true;
            break;
        }
        
        // Check if we're stuck
        if (std::abs(maxForce - prevMaxForce) < 0.01 * convergenceTolerance) {
            stuckCount++;
            if (stuckCount > 3) {
                // Accept current state if close enough
                if (maxForce < convergenceTolerance * 2.0) {
                    converged = true;
                    break;
                }
            }
        } else {
            stuckCount = 0;
        }
        prevMaxForce = maxForce;
        
        // Debug output
        if (iter == maxIterations - 1 && !converged) {
            std::cerr << "FBP: Final iteration " << iter << ", maxForce = " << maxForce 
                      << " kJ/mol/nm (tolerance = " << convergenceTolerance << ")" << std::endl;
        }
        
        // Adaptive damping based on convergence behavior
        if (iter > 0) {
            if (maxForce > convergenceTolerance * 50) {
                // Very large forces - use heavy damping
                dampingFactor = 0.3;
            } else if (maxForce > convergenceTolerance * 10) {
                // Large forces - reduce damping
                dampingFactor = std::min(dampingFactor * 0.9, 0.5);
            } else if (maxForce < convergenceTolerance * 2) {
                // Close to convergence - increase damping for faster convergence
                dampingFactor = std::min(dampingFactor * 1.1, 0.95);
            }
        }
    }
    
    return converged;
}

void DrudeForce::calculateFixedElectricField(const model::MCState& state,
                                             std::vector<Vec3>& fixedField) {
    // Calculate electric field from non-Drude atoms only
    fixedField.assign(particles.size(), Vec3(0.0, 0.0, 0.0));
    
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        const auto& drude = state.atoms[particle.drudeIndex];
        
        // Find which residue this Drude belongs to
        int drudeResIdx = -1;
        for (int resIdx = 0; resIdx < state.activeResidueCount; resIdx++) {
            const auto& res = state.residues[resIdx];
            if (particle.parentIndex >= res.atomStart && 
                particle.parentIndex < res.atomStart + res.atomCount) {
                drudeResIdx = resIdx;
                break;
            }
        }
        
        // Sum field from all non-Drude atoms
        for (int j = 0; j < state.activeAtomCount; j++) {
            const auto& atom = state.atoms[j];
            
            // Skip self
            if (j == particle.drudeIndex) continue;
            
            // Skip other Drude particles
            bool isDrude = false;
            for (const auto& p : particles) {
                if (j == p.drudeIndex) {
                    isDrude = true;
                    break;
                }
            }
            if (isDrude) continue;
            
            // Check intramolecular exclusion
            if (drudeResIdx >= 0) {
                const auto& res = state.residues[drudeResIdx];
                if (j >= res.atomStart && j < res.atomStart + res.atomCount) {
                    continue;  // Same residue - exclude
                }
            }
            
            // Calculate field contribution
            // E = k*q/r^2 * r_hat where r_hat points from source to field point
            double dx = drude.x - atom.x;  // r_field - r_source
            double dy = drude.y - atom.y;
            double dz = drude.z - atom.z;
            
            applyPBC(dx, dy, dz, state.info.box[0], state.info.box[1], state.info.box[2]);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 0.01*0.01) continue;
            
            double r = std::sqrt(r2);
            // Electric field: E = k*q/r^2 * r_hat
            // where r_hat = (r_field - r_source) / |r_field - r_source|
            double factor = 138.935 * atom.charge / (r2 * r);  // k*q/r^3
            
            fixedField[i].x += factor * dx;  // E_x = k*q/r^3 * dx
            fixedField[i].y += factor * dy;
            fixedField[i].z += factor * dz;
        }
    }
}

void DrudeForce::calculateDrudeForces(const model::MCState& state,
                                     std::vector<Vec3>& forces,
                                     const std::vector<Vec3>& fixedField) {
    // Calculate total forces on Drude particles including:
    // 1. Fixed field force
    // 2. Drude-Drude interactions
    // 3. Thole screening corrections
    
    forces = fixedField;  // Start with fixed field
    
    // Add Drude-Drude interactions
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle_i = particles[i];
        const auto& drude_i = state.atoms[particle_i.drudeIndex];
        
        for (size_t j = i + 1; j < particles.size(); j++) {
            const auto& particle_j = particles[j];
            const auto& drude_j = state.atoms[particle_j.drudeIndex];
            
            // Check if in different residues
            bool differentRes = true;
            for (int resIdx = 0; resIdx < state.activeResidueCount; resIdx++) {
                const auto& res = state.residues[resIdx];
                bool i_in_res = (particle_i.parentIndex >= res.atomStart && 
                                particle_i.parentIndex < res.atomStart + res.atomCount);
                bool j_in_res = (particle_j.parentIndex >= res.atomStart && 
                                particle_j.parentIndex < res.atomStart + res.atomCount);
                if (i_in_res && j_in_res) {
                    differentRes = false;
                    break;
                }
            }
            
            if (!differentRes) continue;
            
            // Calculate interaction - field at i due to j
            double dx = drude_i.x - drude_j.x;  // r_i - r_j
            double dy = drude_i.y - drude_j.y;
            double dz = drude_i.z - drude_j.z;
            
            applyPBC(dx, dy, dz, state.info.box[0], state.info.box[1], state.info.box[2]);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 0.01*0.01) continue;
            
            double r = std::sqrt(r2);
            
            // Check for Thole screening
            double thole = 0.0;
            for (const auto& pair : screenedPairs) {
                if ((pair.dipole1 == i && pair.dipole2 == j) ||
                    (pair.dipole1 == j && pair.dipole2 == i)) {
                    thole = pair.thole;
                    break;
                }
            }
            
            // Calculate force with optional Thole screening
            double factor;
            if (thole > 0.0) {
                // Thole exponential screening
                double alpha_i = particle_i.polarizability;
                double alpha_j = particle_j.polarizability;
                double u = r / std::pow(alpha_i * alpha_j, 1.0/6.0);
                double screening = 1.0 - (1.0 + thole*u/2.0) * std::exp(-thole*u);
                factor = 138.935 * particle_j.charge * screening / (r2 * r);
            } else {
                factor = 138.935 * particle_j.charge / (r2 * r);
            }
            
            // Update forces (Newton's third law)
            forces[i].x += factor * dx;
            forces[i].y += factor * dy;
            forces[i].z += factor * dz;
            
            // Equal and opposite force on j (scaled by charge ratio)
            forces[j].x -= factor * dx;
            forces[j].y -= factor * dy;
            forces[j].z -= factor * dz;
        }
    }
    
    // Convert electric field to force by multiplying by charge
    for (size_t i = 0; i < particles.size(); i++) {
        forces[i] = forces[i] * particles[i].charge;
    }
}

double DrudeForce::calculateEnergyFBP(model::MCState& state) {
    // Use Force Balance Predictor to optimize Drude positions
    bool converged = minimizeDrudePositionsWithFBP(state);
    
    if (!converged) {
        std::cerr << "Warning: FBP did not converge" << std::endl;
    }
    
    // Calculate final energy
    return calculateEnergyDirect(state);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc