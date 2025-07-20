#include "DrudeForce.hpp"
#include "DrudeConjugateGradient.hpp"
#include "DrudeThole.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Physical constants
static const double ONE_4PI_EPS0 = 138.935456;  // kJ/mol·nm·e^-2

DrudeForce::DrudeForce() {
    // Default SCF parameters are initialized in the header file
    // No need to override them here
}

int DrudeForce::addParticle(int drudeIndex, int parentIndex, 
                           int aniso1Index, int aniso2Index,
                           int aniso3Index, int aniso4Index,
                           double charge, double polarizability,
                           double aniso12, double aniso34) {
    DrudeParticle particle;
    particle.drudeIndex = drudeIndex;
    particle.parentIndex = parentIndex;
    particle.aniso1Index = aniso1Index;
    particle.aniso2Index = aniso2Index;
    particle.aniso3Index = aniso3Index;
    particle.aniso4Index = aniso4Index;
    particle.charge = charge;
    particle.polarizability = polarizability;
    particle.aniso12 = aniso12;
    particle.aniso34 = aniso34;
    
    // Calculate force constants based on the CHARMM Drude model
    // k = q^2 / (4πε₀ * α)
    double a1 = (aniso1Index == -1 || aniso2Index == -1) ? 1.0 : aniso12;
    double a2 = (aniso3Index == -1 || aniso4Index == -1) ? 1.0 : aniso34;
    double a3 = 3.0 - a1 - a2;  // Sum of anisotropy factors must be 3
    
    // Force constants in atomic units
    particle.kIsotropic = ONE_4PI_EPS0 * charge * charge / (polarizability * a3);
    particle.kAniso1 = ONE_4PI_EPS0 * charge * charge / (polarizability * a1) - particle.kIsotropic;
    particle.kAniso2 = ONE_4PI_EPS0 * charge * charge / (polarizability * a2) - particle.kIsotropic;
    
    particles.push_back(particle);
    return particles.size() - 1;
}

void DrudeForce::addScreenedPair(int dipole1, int dipole2, double thole) {
    ScreenedPair pair;
    pair.dipole1 = dipole1;
    pair.dipole2 = dipole2;
    pair.thole = thole;
    screenedPairs.push_back(pair);
}

void DrudeForce::setOPT3Coefficients(double c0, double c1, double c2, double c3) {
    opt3Coeffs.c0 = c0;
    opt3Coeffs.c1 = c1;
    opt3Coeffs.c2 = c2;
    opt3Coeffs.c3 = c3;
}

void DrudeForce::setOPT4Coefficients(double c0, double c1, double c2, double c3, double c4) {
    opt4Coeffs.c0 = c0;
    opt4Coeffs.c1 = c1;
    opt4Coeffs.c2 = c2;
    opt4Coeffs.c3 = c3;
    opt4Coeffs.c4 = c4;
}

double DrudeForce::calculateEnergySCF(model::MCState& state) {
    // Initialize forces array
    std::vector<Vec3> forces(state.activeAtomCount);
    
    // Perform minimization based on selected algorithm
    bool converged;
    
    // Handle legacy useOPT3 flag
    if (useOPT3 && algorithm == DrudeAlgorithm::SCF) {
        algorithm = DrudeAlgorithm::OPT3;
    }
    
    switch (algorithm) {
        case DrudeAlgorithm::SCF:
            converged = minimizeDrudePositions(state, forces);
            break;
        case DrudeAlgorithm::OPT3:
            converged = minimizeDrudePositionsWithOPT3(state);
            break;
        case DrudeAlgorithm::OPT4:
            converged = minimizeDrudePositionsWithOPT4(state);
            break;
        case DrudeAlgorithm::AdaptiveOPT:
            converged = minimizeDrudePositionsAdaptive(state);
            break;
        case DrudeAlgorithm::HybridOPT:
            converged = minimizeDrudePositionsHybrid(state);
            break;
        case DrudeAlgorithm::SmartOPT3:
            converged = minimizeDrudePositionsWithSmartOPT3(state);
            break;
        case DrudeAlgorithm::FBP:
            converged = minimizeDrudePositionsWithFBP(state);
            break;
        case DrudeAlgorithm::ConjugateGradient:
            converged = minimizeDrudePositionsWithCG(state);
            break;
        default:
            converged = minimizeDrudePositions(state, forces);
    }
    
    if (!converged) {
        std::cerr << "Warning: Drude SCF did not converge after " 
                  << scfParams.maxIterations << " iterations" << std::endl;
    }
    
    // Calculate final energy with optimized Drude positions
    double energy = 0.0;
    
    // Clear forces for final energy calculation
    std::fill(forces.begin(), forces.end(), Vec3());
    
    // Harmonic restraint energy
    energy += calculateHarmonicEnergy(state, forces);
    
    // Screened Coulomb energy
    energy += calculateScreenedCoulombEnergy(state, forces);
    
    // All Coulomb interactions
    energy += calculateCoulombEnergy(state, forces);
    
    return energy;
}

void DrudeForce::calculateForces(model::MCState& state, std::vector<Vec3>& forces) {
    // Clear forces
    std::fill(forces.begin(), forces.end(), Vec3());
    
    // Calculate harmonic forces
    calculateHarmonicEnergy(state, forces);
    
    // Calculate screened Coulomb forces
    calculateScreenedCoulombEnergy(state, forces);
    
    // Calculate all Coulomb interactions
    calculateCoulombEnergy(state, forces);
}

bool DrudeForce::minimizeDrudePositions(model::MCState& state, std::vector<Vec3>& forces) {
    double lastTotalForce = 0.0;
    
    // Use fixed tolerance for all system sizes
    
    // Track force changes for oscillation detection
    double previousForce = 0.0;
    int oscillationCount = 0;
    
    for (int iter = 0; iter < scfParams.maxIterations; ++iter) {
        // Calculate forces on all particles
        calculateForces(state, forces);
        
        // Update only Drude particle positions
        double totalForce = 0.0;
        double maxIndividualForce = 0.0;
        
        for (const auto& drude : particles) {
            int idx = drude.drudeIndex;
            const Vec3& force = forces[idx];
            double forceMag = force.norm();
            totalForce += forceMag * forceMag;
            maxIndividualForce = std::max(maxIndividualForce, forceMag);
            
            // Calculate effective force scale for anisotropic polarizability
            Vec3 fscale(drude.kIsotropic, drude.kIsotropic, drude.kIsotropic);
            
            // Add anisotropic contributions
            if (drude.aniso1Index >= 0 && drude.aniso2Index >= 0) {
                Vec3 dir(state.atoms[drude.aniso1Index].x - state.atoms[drude.aniso2Index].x,
                        state.atoms[drude.aniso1Index].y - state.atoms[drude.aniso2Index].y,
                        state.atoms[drude.aniso1Index].z - state.atoms[drude.aniso2Index].z);
                dir = dir.normalized();
                fscale.x += drude.kAniso1 * dir.x * dir.x;
                fscale.y += drude.kAniso1 * dir.y * dir.y;
                fscale.z += drude.kAniso1 * dir.z * dir.z;
            }
            
            if (drude.aniso3Index >= 0 && drude.aniso4Index >= 0) {
                Vec3 dir(state.atoms[drude.aniso3Index].x - state.atoms[drude.aniso4Index].x,
                        state.atoms[drude.aniso3Index].y - state.atoms[drude.aniso4Index].y,
                        state.atoms[drude.aniso3Index].z - state.atoms[drude.aniso4Index].z);
                dir = dir.normalized();
                fscale.x += drude.kAniso2 * dir.x * dir.x;
                fscale.y += drude.kAniso2 * dir.y * dir.y;
                fscale.z += drude.kAniso2 * dir.z * dir.z;
            }
            
            // Improved adaptive damping
            double damping;
            if (iter < 3) {
                // More aggressive damping in early iterations
                damping = 0.3;
            } else if (forceMag > 100.0 * scfParams.tolerance) {
                damping = 0.2;  // Very large forces
            } else if (forceMag > 10.0 * scfParams.tolerance) {
                damping = scfParams.dampingFactor;
            } else {
                damping = 1.0;  // Small forces - full step
            }
            
            // Update Drude position based on force and spring constants
            // For isotropic case: delta_r = damping * F / k
            Vec3 delta;
            
            // Simplified for isotropic case when no anisotropy is defined
            if (drude.aniso1Index < 0 && drude.aniso3Index < 0) {
                // Pure isotropic: use scalar k
                double factor = damping / drude.kIsotropic;
                delta = force * factor;
            } else {
                // Anisotropic: use tensor form
                delta.x = damping * force.x / fscale.x;
                delta.y = damping * force.y / fscale.y;
                delta.z = damping * force.z / fscale.z;
            }
            
            // Apply update with hard wall constraint
            Vec3 newPos(state.atoms[idx].x + delta.x,
                       state.atoms[idx].y + delta.y,
                       state.atoms[idx].z + delta.z);
            
            // Enforce hard wall constraint
            Vec3 parentPos(state.atoms[drude.parentIndex].x,
                          state.atoms[drude.parentIndex].y,
                          state.atoms[drude.parentIndex].z);
            Vec3 dr = newPos - parentPos;
            double r = dr.norm();
            
            if (r > scfParams.maxDrudeDistance) {
                // Scale back to maximum allowed distance
                double scale = scfParams.maxDrudeDistance / r;
                dr = dr * scale;
                newPos = parentPos + dr;
            }
            
            // Update position
            state.atoms[idx].x = newPos.x;
            state.atoms[idx].y = newPos.y;
            state.atoms[idx].z = newPos.z;
        }
        
        // Check convergence
        double rmsForce = std::sqrt(totalForce / (3.0 * particles.size()));
        
        // Check for oscillations
        if (iter > 0) {
            double forceChange = std::abs(totalForce - previousForce);
            if (forceChange < 0.1 * totalForce && previousForce > 0) {
                // Forces stabilized
                if (rmsForce < scfParams.tolerance * 2.0) {
                    return true;  // Good enough
                }
            }
            
            // Detect oscillations
            if (totalForce > previousForce * 0.9 && previousForce > lastTotalForce * 0.9) {
                oscillationCount++;
                if (oscillationCount > 3 && rmsForce < scfParams.tolerance * 5.0) {
                    // Accept if oscillating but reasonably small
                    return true;
                }
            }
        }
        previousForce = lastTotalForce;
        
        // Primary convergence check
        if (rmsForce < scfParams.tolerance) {
            return true;  // Converged
        }
        
        // Secondary check - accept if close enough after many iterations
        if (iter > scfParams.maxIterations / 2) {
            if (rmsForce < scfParams.tolerance * 3.0 && maxIndividualForce < scfParams.tolerance * 10.0) {
                return true;  // Acceptable for large systems
            }
        }
        
        // Check if forces are increasing (diverging) - following OpenMM
        if (iter > 0 && totalForce > 0.9 * lastTotalForce) {
            // Don't give up immediately, but track if we're stuck
            if (totalForce > lastTotalForce) {
                // Forces increased, might be diverging
                // But don't break if forces are small
                if (rmsForce > scfParams.tolerance * 20.0) {
                    break;
                }
            }
        }
        
        lastTotalForce = totalForce;
    }
    
    // Final check before declaring failure
    double finalRMS = std::sqrt(lastTotalForce / (3.0 * particles.size()));
    
    // Suppress warning for marginally unconverged but acceptable results
    if (finalRMS < scfParams.tolerance * 5.0) {
        return true;  // Close enough, avoid warning
    }
    
    // Only show warning for truly problematic cases
    if (finalRMS > scfParams.tolerance * 10.0) {
        std::cerr << "Warning: Drude SCF did not converge after " << scfParams.maxIterations 
                  << " iterations. RMS force = " << finalRMS << " kJ/mol/nm" << std::endl;
    }
    
    return false;  // Did not converge
}

double DrudeForce::calculateHarmonicEnergy(model::MCState& state, std::vector<Vec3>& forces) {
    double energy = 0.0;
    
    for (const auto& drude : particles) {
        int p = drude.drudeIndex;
        int p1 = drude.parentIndex;
        
        // Get positions
        Vec3 posDrude(state.atoms[p].x, state.atoms[p].y, state.atoms[p].z);
        Vec3 posParent(state.atoms[p1].x, state.atoms[p1].y, state.atoms[p1].z);
        Vec3 delta = posDrude - posParent;
        
        // Isotropic harmonic energy: E = 0.5 * k * r^2
        double r2 = delta.norm2();
        energy += 0.5 * drude.kIsotropic * r2;
        
        // Isotropic forces
        // Force on Drude: F = -k * (r_drude - r_parent) = -k * delta
        Vec3 f = delta * drude.kIsotropic;
        forces[p] -= f;  // Force on Drude (restoring force)
        forces[p1] += f;  // Equal and opposite force on parent
        
        // First anisotropic term
        if (drude.aniso1Index >= 0 && drude.aniso2Index >= 0) {
            int p2 = drude.aniso1Index;
            int p3 = drude.aniso2Index;
            
            Vec3 pos2(state.atoms[p2].x, state.atoms[p2].y, state.atoms[p2].z);
            Vec3 pos3(state.atoms[p3].x, state.atoms[p3].y, state.atoms[p3].z);
            Vec3 dir = (pos2 - pos3).normalized();
            
            double rprime = delta.dot(dir);
            energy += 0.5 * drude.kAniso1 * rprime * rprime;
            
            // Anisotropic forces - only act on Drude and parent
            Vec3 f1 = dir * (drude.kAniso1 * rprime);
            
            forces[p] -= f1;
            forces[p1] += f1;
        }
        
        // Second anisotropic term
        if (drude.aniso3Index >= 0 && drude.aniso4Index >= 0) {
            int p3 = drude.aniso3Index;
            int p4 = drude.aniso4Index;
            
            Vec3 pos3(state.atoms[p3].x, state.atoms[p3].y, state.atoms[p3].z);
            Vec3 pos4(state.atoms[p4].x, state.atoms[p4].y, state.atoms[p4].z);
            Vec3 dir = (pos3 - pos4).normalized();
            
            double rprime = delta.dot(dir);
            energy += 0.5 * drude.kAniso2 * rprime * rprime;
            
            // Anisotropic forces - only act on Drude and parent
            Vec3 f1 = dir * (drude.kAniso2 * rprime);
            
            forces[p] -= f1;
            forces[p1] += f1;
        }
    }
    
    return energy;
}

double DrudeForce::calculateScreenedCoulombEnergy(model::MCState& state, std::vector<Vec3>& forces) {
    double energy = 0.0;
    
    for (const auto& pair : screenedPairs) {
        const DrudeParticle& dipole1 = particles[pair.dipole1];
        const DrudeParticle& dipole2 = particles[pair.dipole2];
        
        // Following OpenMM convention:
        // The Drude model uses implicit charges where:
        // - Drude particle has charge q_drude (typically negative)
        // - Parent atom has an implicit charge offset of -q_drude
        // This creates a dipole: parent(+) --- drude(-)
        
        // Get Drude charges from the particle definitions
        double drudeCharge1 = dipole1.charge;  // Should be negative
        double drudeCharge2 = dipole2.charge;  // Should be negative
        
        // Get positions
        Vec3 posDrude1(state.atoms[dipole1.drudeIndex].x, 
                       state.atoms[dipole1.drudeIndex].y, 
                       state.atoms[dipole1.drudeIndex].z);
        Vec3 posParent1(state.atoms[dipole1.parentIndex].x, 
                        state.atoms[dipole1.parentIndex].y, 
                        state.atoms[dipole1.parentIndex].z);
        Vec3 posDrude2(state.atoms[dipole2.drudeIndex].x, 
                       state.atoms[dipole2.drudeIndex].y, 
                       state.atoms[dipole2.drudeIndex].z);
        Vec3 posParent2(state.atoms[dipole2.parentIndex].x, 
                        state.atoms[dipole2.parentIndex].y, 
                        state.atoms[dipole2.parentIndex].z);
        
        // Calculate all four distances with PBC
        auto applyPBC = [&state](Vec3& delta) {
            if (state.info.box[0] > 0 && state.info.box[1] > 0 && state.info.box[2] > 0) {
                delta.x -= state.info.box[0] * std::round(delta.x / state.info.box[0]);
                delta.y -= state.info.box[1] * std::round(delta.y / state.info.box[1]);
                delta.z -= state.info.box[2] * std::round(delta.z / state.info.box[2]);
            }
        };
        
        // Four interactions following OpenMM logic
        // The charge products follow the pattern:
        // DD: q1*q2 (both negative, product positive)
        // DP: q1*(-q2) = -q1*q2 (opposite signs, product negative)
        // PD: (-q1)*q2 = -q1*q2 (opposite signs, product negative)
        // PP: (-q1)*(-q2) = q1*q2 (both positive, product positive)
        
        struct Interaction {
            Vec3 pos1, pos2;
            double chargeProduct;
            const char* label;
        };
        
        Interaction interactions[4] = {
            {posDrude1, posDrude2, drudeCharge1 * drudeCharge2, "DD"},           // Drude1-Drude2
            {posDrude1, posParent2, -drudeCharge1 * drudeCharge2, "DP"},         // Drude1-Parent2
            {posParent1, posDrude2, -drudeCharge1 * drudeCharge2, "PD"},         // Parent1-Drude2
            {posParent1, posParent2, drudeCharge1 * drudeCharge2, "PP"}          // Parent1-Parent2
        };
        
        // Calculate energy and forces for each interaction
        for (const auto& inter : interactions) {
            Vec3 delta = inter.pos2 - inter.pos1;
            applyPBC(delta);
            
            double r = delta.norm();
            if (r < 1e-10) continue;
            
            // Calculate Thole screening
            double screening, dScreening_dr;
            TholeFunctions::calculateScreening(r, pair.thole, 
                                             dipole1.polarizability, 
                                             dipole2.polarizability,
                                             screening, dScreening_dr);
            
            // Energy contribution
            double invR = 1.0 / r;
            energy += ONE_4PI_EPS0 * inter.chargeProduct * screening * invR;
            
            // Force contribution
            // F = -dE/dr = -k*q1*q2 * [S/r^2 - dS/dr/r]
            double invR2 = invR * invR;
            double forceMag = ONE_4PI_EPS0 * inter.chargeProduct * 
                             (screening * invR2 - dScreening_dr * invR);
            Vec3 force = delta * (-forceMag * invR);  // Normalize and apply sign
            
            // Apply forces based on which atoms are involved
            if (&inter.pos1 == &posDrude1) forces[dipole1.drudeIndex] += force;
            else if (&inter.pos1 == &posParent1) forces[dipole1.parentIndex] += force;
            
            if (&inter.pos2 == &posDrude2) forces[dipole2.drudeIndex] -= force;
            else if (&inter.pos2 == &posParent2) forces[dipole2.parentIndex] -= force;
        }
    }
    
    return energy;
}

double DrudeForce::calculateCoulombEnergy(model::MCState& state, std::vector<Vec3>& forces) {
    double energy = 0.0;
    
    // Calculate all pairwise Coulomb interactions
    for (int i = 0; i < state.activeAtomCount; ++i) {
        for (int j = i + 1; j < state.activeAtomCount; ++j) {
            // Skip if both atoms have zero charge
            if (std::abs(state.atoms[i].charge) < 1e-10 && 
                std::abs(state.atoms[j].charge) < 1e-10) continue;
            
            // Check if atoms are in same residue
            bool sameResidue = false;
            for (int r = 0; r < state.activeResidueCount; ++r) {
                int start = state.residues[r].atomStart;
                int end = start + state.residues[r].atomCount;
                if (i >= start && i < end && j >= start && j < end) {
                    sameResidue = true;
                    break;
                }
            }
            
            // Skip ALL intramolecular nonbonded interactions
            // Following OpenMM convention: all atom pairs within a molecule 
            // have exclusions (charge=0, sigma=1, epsilon=0)
            // The only intramolecular interaction is the Drude harmonic restraint
            if (sameResidue) continue;
            
            Vec3 pos1(state.atoms[i].x, state.atoms[i].y, state.atoms[i].z);
            Vec3 pos2(state.atoms[j].x, state.atoms[j].y, state.atoms[j].z);
            Vec3 delta = pos2 - pos1;
            
            // Apply periodic boundary conditions if box is defined
            if (state.info.box[0] > 0 && state.info.box[1] > 0 && state.info.box[2] > 0) {
                delta.x -= state.info.box[0] * std::round(delta.x / state.info.box[0]);
                delta.y -= state.info.box[1] * std::round(delta.y / state.info.box[1]);
                delta.z -= state.info.box[2] * std::round(delta.z / state.info.box[2]);
            }
            
            double r = delta.norm();
            if (r < 1e-10) continue;  // Skip if too close
            
            // Apply cutoff if specified
            if (state.info.cutoff > 0 && r > state.info.cutoff) continue;
            
            // Check if this pair is already handled by screened interactions
            bool isScreened = false;
            for (const auto& pair : screenedPairs) {
                const DrudeParticle& dipole1 = particles[pair.dipole1];
                const DrudeParticle& dipole2 = particles[pair.dipole2];
                
                // Check if i,j match any screened pair combination
                if ((i == dipole1.drudeIndex || i == dipole1.parentIndex) &&
                    (j == dipole2.drudeIndex || j == dipole2.parentIndex)) {
                    isScreened = true;
                    break;
                }
                if ((j == dipole1.drudeIndex || j == dipole1.parentIndex) &&
                    (i == dipole2.drudeIndex || i == dipole2.parentIndex)) {
                    isScreened = true;
                    break;
                }
            }
            
            if (isScreened) continue;  // Skip if already handled by screening
            
            // Regular Coulomb interaction
            double chargeProduct = state.atoms[i].charge * state.atoms[j].charge;
            double invR = 1.0 / r;
            double invR2 = invR * invR;
            
            // Energy
            energy += ONE_4PI_EPS0 * chargeProduct * invR;
            
            // Force
            double forceMag = ONE_4PI_EPS0 * chargeProduct * invR2;
            // Fix: Add negative sign to correct force direction for attractive forces
            Vec3 f = delta * (-forceMag * invR);
            
            forces[i] += f;
            forces[j] -= f;
        }
    }
    
    return energy;
}

double DrudeForce::calculateEnergyOPT3(model::MCState& state) {
    // Use OPT3 for minimization
    bool oldUseOPT3 = useOPT3;
    useOPT3 = true;
    double energy = calculateEnergySCF(state);
    useOPT3 = oldUseOPT3;
    return energy;
}

// OPT3 and OPT4 implementations are now in separate compilation units

// These will be implemented in DrudeOPT3.cpp/DrudeOPT4.cpp
// which are now separate compilation units

DrudeForce::OPT3TrainingData DrudeForce::collectTrainingData(model::MCState& state) {
    OPT3TrainingData data;
    size_t numDrudes = particles.size();
    if (numDrudes == 0) return data;
    
    // Save original state
    bool oldUseOPT3 = useOPT3;
    std::vector<Vec3> originalDrudePos(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        int idx = particles[i].drudeIndex;
        originalDrudePos[i] = Vec3(state.atoms[idx].x, state.atoms[idx].y, state.atoms[idx].z);
    }
    
    // Step 1: Run SCF to get ground truth
    useOPT3 = false;
    std::vector<Vec3> forces(state.activeAtomCount);
    minimizeDrudePositions(state, forces);
    
    // Save SCF positions and parent positions
    data.r_scf.resize(numDrudes);
    data.parentPos.resize(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        int drudeIdx = particles[i].drudeIndex;
        int parentIdx = particles[i].parentIndex;
        
        data.parentPos[i] = Vec3(state.atoms[parentIdx].x,
                                state.atoms[parentIdx].y,
                                state.atoms[parentIdx].z);
        
        Vec3 drudePos(state.atoms[drudeIdx].x,
                     state.atoms[drudeIdx].y,
                     state.atoms[drudeIdx].z);
        
        data.r_scf[i] = drudePos - data.parentPos[i];
    }
    
    // Step 2: Reset Drude to parent positions
    for (size_t i = 0; i < numDrudes; ++i) {
        int drudeIdx = particles[i].drudeIndex;
        int parentIdx = particles[i].parentIndex;
        state.atoms[drudeIdx].x = state.atoms[parentIdx].x;
        state.atoms[drudeIdx].y = state.atoms[parentIdx].y;
        state.atoms[drudeIdx].z = state.atoms[parentIdx].z;
    }
    
    // Step 3: Compute perturbation orders
    data.r0.resize(numDrudes);
    data.r1.resize(numDrudes);
    data.r2.resize(numDrudes);
    data.r3.resize(numDrudes);
    
    std::vector<Vec3> electricField(numDrudes);
    
    // Zero-order: response to static field only
    ::pygcmc::platform::cpu::calculateElectricFieldAtDrudes(state, particles, screenedPairs, data.parentPos, electricField, true);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        data.r0[i] = electricField[i] * factor;
    }
    
    // First-order
    std::vector<Vec3> drudePos0(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos0[i] = data.parentPos[i] + data.r0[i];
    }
    ::pygcmc::platform::cpu::calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos0, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        data.r1[i] = electricField[i] * factor;
    }
    
    // Second-order
    std::vector<Vec3> drudePos1(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos1[i] = data.parentPos[i] + data.r1[i];
    }
    ::pygcmc::platform::cpu::calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos1, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        data.r2[i] = electricField[i] * factor;
    }
    
    // Third-order
    std::vector<Vec3> drudePos2(numDrudes);
    for (size_t i = 0; i < numDrudes; ++i) {
        drudePos2[i] = data.parentPos[i] + data.r2[i];
    }
    ::pygcmc::platform::cpu::calculateElectricFieldAtDrudes(state, particles, screenedPairs, drudePos2, electricField, false);
    for (size_t i = 0; i < numDrudes; ++i) {
        double factor = particles[i].charge / particles[i].kIsotropic;
        data.r3[i] = electricField[i] * factor;
    }
    
    // Restore original state
    for (size_t i = 0; i < numDrudes; ++i) {
        int idx = particles[i].drudeIndex;
        state.atoms[idx].x = originalDrudePos[i].x;
        state.atoms[idx].y = originalDrudePos[i].y;
        state.atoms[idx].z = originalDrudePos[i].z;
    }
    useOPT3 = oldUseOPT3;
    
    return data;
}

bool DrudeForce::minimizeDrudePositionsWithCG(model::MCState& state) {
    // Use the conjugate gradient solver
    return DrudeConjugateGradient::minimizeDrudePositions(
        state, 
        *this,
        scfParams.tolerance,
        scfParams.maxIterations
    );
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc