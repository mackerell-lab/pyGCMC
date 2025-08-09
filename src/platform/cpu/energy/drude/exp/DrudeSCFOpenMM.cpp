#include "DrudeSCFOpenMM.hpp"
#include "../DrudeStructures.hpp"  // Contains DrudeConstants
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

void DrudeSCFOpenMM::applyPBC(double& dx, double& dy, double& dz, 
                               const std::array<double, 3>& box) {
    if (box[0] > 0) dx -= box[0] * std::round(dx / box[0]);
    if (box[1] > 0) dy -= box[1] * std::round(dy / box[1]);
    if (box[2] > 0) dz -= box[2] * std::round(dz / box[2]);
}

double DrudeSCFOpenMM::tholeS3(double r, double alpha_i, double alpha_j, double thole) const {
    if (thole <= 1e-14 || alpha_i <= 1e-14 || alpha_j <= 1e-14) {
        return 1.0;  // No screening
    }
    
    double alpha_eff = std::pow(alpha_i * alpha_j, 1.0/6.0);
    double u = thole * r / alpha_eff;
    
    if (u > 50.0) {
        return 1.0;  // Screening is negligible for large u
    }
    
    double exp_u = std::exp(-u);
    // Thole S3 function for dipole field (1/r^3 kernel)
    // S3(u) = 1 - exp(-u) * (1 + u + 0.5*u^2)
    return 1.0 - exp_u * (1.0 + u + 0.5 * u * u);
}

bool DrudeSCFOpenMM::optimize(model::MCState& state,
                              const std::vector<DrudeParticle>& particles,
                              const std::vector<ScreenedPair>& pairs,
                              const DrudeSCFParams& params) {
    if (particles.empty()) {
        return true;
    }
    
    // Initialize electric field vectors
    std::vector<Vec3> electricField(particles.size(), {0.0, 0.0, 0.0});
    
    // Best-so-far tracking
    double bestEnergy = 1e300;
    std::vector<DrudePosition> bestPositions(particles.size());
    
    // Adaptive damping
    double damping = params.dampingFactor;
    double previousEnergy = 1e300;
    int stallCount = 0;
    
    // SCF iteration
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        // Clear electric field
        std::fill(electricField.begin(), electricField.end(), Vec3{0.0, 0.0, 0.0});
        
        // Calculate total electric field
        calculateElectricField(state, particles, pairs, electricField);
        
        // Check convergence based on forces
        double maxForce = 0.0;
        for (size_t i = 0; i < particles.size(); ++i) {
            const auto& p = particles[i];
            if (p.polarizability < 1e-14) continue;  // Skip frozen particles
            
            const auto& drude = state.atoms[p.drudeIndex];
            const auto& parent = state.atoms[p.parentIndex];
            
            double dx = drude.x - parent.x;
            double dy = drude.y - parent.y;
            double dz = drude.z - parent.z;
            
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(dx, dy, dz, box);
            
            // Force = q*E - k*d
            double fx = p.charge * electricField[i][0] - p.kSpring * dx;
            double fy = p.charge * electricField[i][1] - p.kSpring * dy;
            double fz = p.charge * electricField[i][2] - p.kSpring * dz;
            
            double force2 = fx*fx + fy*fy + fz*fz;
            maxForce = std::max(maxForce, std::sqrt(force2));
        }
        
        // Check convergence
        if (maxForce < params.tolerance) {
            return true;  // Converged
        }
        
        // Update positions
        double maxStep = 0.002;  // 0.002 nm max step per iteration
        double hardWall = params.enableHardWall ? params.maxDrudeDistance : 1e9;
        updateDrudePositions(state, particles, electricField, damping, maxStep, hardWall);
        
        // Calculate current energy
        double currentEnergy = calculateSpringEnergy(state, particles);
        
        // Track best-so-far
        if (currentEnergy < bestEnergy) {
            bestEnergy = currentEnergy;
            for (size_t i = 0; i < particles.size(); ++i) {
                const auto& p = particles[i];
                const auto& drude = state.atoms[p.drudeIndex];
                const auto& parent = state.atoms[p.parentIndex];
                bestPositions[i].dx = drude.x - parent.x;
                bestPositions[i].dy = drude.y - parent.y;
                bestPositions[i].dz = drude.z - parent.z;
            }
            stallCount = 0;
        } else {
            stallCount++;
        }
        
        // Adaptive damping
        if (currentEnergy < previousEnergy) {
            damping *= 0.95;  // Reduce damping if improving
            damping = std::max(damping, 0.05);
        } else if (stallCount > 5) {
            damping *= 1.1;  // Increase damping if stalled
            damping = std::min(damping, 0.5);
            stallCount = 0;
        }
        
        previousEnergy = currentEnergy;
    }
    
    // Revert to best-so-far positions
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        auto& drude = state.atoms[p.drudeIndex];
        const auto& parent = state.atoms[p.parentIndex];
        drude.x = parent.x + bestPositions[i].dx;
        drude.y = parent.y + bestPositions[i].dy;
        drude.z = parent.z + bestPositions[i].dz;
    }
    
    return false;  // Did not converge but returned best found
}

void DrudeSCFOpenMM::calculateElectricField(const model::MCState& state,
                                            const std::vector<DrudeParticle>& particles,
                                            const std::vector<ScreenedPair>& pairs,
                                            std::vector<Vec3>& electricField) const {
    calculateExternalField(state, particles, pairs, electricField);
    calculateInducedField(state, particles, pairs, electricField);
}

void DrudeSCFOpenMM::calculateExternalField(const model::MCState& state,
                                            const std::vector<DrudeParticle>& particles,
                                            const std::vector<ScreenedPair>& pairs,
                                            std::vector<Vec3>& electricField) const {
    // Build exclusion list: for each Drude, list of parent atoms from screened pairs
    std::vector<std::vector<int>> excludeParents(particles.size());
    for (const auto& pair : pairs) {
        if (static_cast<size_t>(pair.dipole1) < particles.size() && 
            static_cast<size_t>(pair.dipole2) < particles.size()) {
            excludeParents[pair.dipole1].push_back(particles[pair.dipole2].parentIndex);
            excludeParents[pair.dipole2].push_back(particles[pair.dipole1].parentIndex);
        }
    }
    
    // Calculate external field from non-Drude charges
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& drudeAtom = state.atoms[particle.drudeIndex];
        
        Vec3 field = {0.0, 0.0, 0.0};
        
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self and parent
            if (j == particle.drudeIndex || j == particle.parentIndex) continue;
            
            // Skip all other Drude particles
            bool isDrude = false;
            for (const auto& otherParticle : particles) {
                if (j == otherParticle.drudeIndex) {
                    isDrude = true;
                    break;
                }
            }
            if (isDrude) continue;
            
            // Skip partner parents (to avoid double counting with induced field)
            if (std::find(excludeParents[i].begin(), excludeParents[i].end(), j) 
                != excludeParents[i].end()) {
                continue;
            }
            
            // Skip intramolecular interactions if needed
            if (inSameMolecule(particle.drudeIndex, j, state)) {
                // For same molecule, might want to skip certain atoms
                // For now, skip if in same molecule
                continue;
            }
            
            const auto& atom = state.atoms[j];
            if (std::abs(atom.charge) < 1e-10) continue;
            
            // Calculate distance
            double dx = drudeAtom.x - atom.x;
            double dy = drudeAtom.y - atom.y;
            double dz = drudeAtom.z - atom.z;
            
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;  // Skip if too close
            
            double r = std::sqrt(r2);
            double r3 = r2 * r;
            
            // Electric field: E = k * q / r^2 * r_hat
            double factor = DrudeConstants::ONE_4PI_EPS0 * atom.charge / r3;
            
            field[0] += factor * dx;
            field[1] += factor * dy;
            field[2] += factor * dz;
        }
        
        electricField[i][0] += field[0];
        electricField[i][1] += field[1];
        electricField[i][2] += field[2];
    }
}

void DrudeSCFOpenMM::calculateInducedField(const model::MCState& state,
                                           const std::vector<DrudeParticle>& particles,
                                           const std::vector<ScreenedPair>& pairs,
                                           std::vector<Vec3>& electricField) const {
    // Calculate field from induced dipoles using dipole tensor with Thole screening
    for (const auto& pair : pairs) {
        if (static_cast<size_t>(pair.dipole1) >= particles.size() || 
            static_cast<size_t>(pair.dipole2) >= particles.size()) {
            continue;
        }
        
        const auto& particle1 = particles[pair.dipole1];
        const auto& particle2 = particles[pair.dipole2];
        
        const auto& parent1 = state.atoms[particle1.parentIndex];
        const auto& drude1 = state.atoms[particle1.drudeIndex];
        const auto& parent2 = state.atoms[particle2.parentIndex];
        const auto& drude2 = state.atoms[particle2.drudeIndex];
        
        // Calculate dipole moments: μ = q_D * (D - P)
        double dx1 = drude1.x - parent1.x;
        double dy1 = drude1.y - parent1.y;
        double dz1 = drude1.z - parent1.z;
        
        double dx2 = drude2.x - parent2.x;
        double dy2 = drude2.y - parent2.y;
        double dz2 = drude2.z - parent2.z;
        
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx1, dy1, dz1, box);
        applyPBC(dx2, dy2, dz2, box);
        
        // Dipole moments (with sign from charge)
        Vec3 mu1 = {particle1.charge * dx1, particle1.charge * dy1, particle1.charge * dz1};
        Vec3 mu2 = {particle2.charge * dx2, particle2.charge * dy2, particle2.charge * dz2};
        
        // Vector from parent1 to parent2
        double rx = parent2.x - parent1.x;
        double ry = parent2.y - parent1.y;
        double rz = parent2.z - parent1.z;
        applyPBC(rx, ry, rz, box);
        
        double r2 = rx*rx + ry*ry + rz*rz;
        if (r2 < 1e-12) continue;
        
        double r = std::sqrt(r2);
        double invr3 = 1.0 / (r2 * r);
        
        // Unit vector from parent1 to parent2
        double nx = rx / r;
        double ny = ry / r;
        double nz = rz / r;
        
        // Calculate Thole screening S3
        double s = tholeS3(r, particle1.polarizability, particle2.polarizability, pair.thole);
        
        // Prefactor with screening
        double prefactor = DrudeConstants::ONE_4PI_EPS0 * s * invr3;
        
        // Dipole tensor field calculation
        // E1 from μ2: E = k * S3/r^3 * [3(μ2·n)n - μ2]
        double dot2 = mu2[0]*nx + mu2[1]*ny + mu2[2]*nz;
        Vec3 E1 = {
            prefactor * (3.0 * dot2 * nx - mu2[0]),
            prefactor * (3.0 * dot2 * ny - mu2[1]),
            prefactor * (3.0 * dot2 * nz - mu2[2])
        };
        
        // E2 from μ1: E = k * S3/r^3 * [3(μ1·n)n - μ1]
        // Note: n points from 1 to 2, so for field at 2 from 1, we use -n
        double dot1 = mu1[0]*(-nx) + mu1[1]*(-ny) + mu1[2]*(-nz);
        Vec3 E2 = {
            prefactor * (3.0 * dot1 * (-nx) - mu1[0]),
            prefactor * (3.0 * dot1 * (-ny) - mu1[1]),
            prefactor * (3.0 * dot1 * (-nz) - mu1[2])
        };
        
        // Add to electric field
        electricField[pair.dipole1][0] += E1[0];
        electricField[pair.dipole1][1] += E1[1];
        electricField[pair.dipole1][2] += E1[2];
        
        electricField[pair.dipole2][0] += E2[0];
        electricField[pair.dipole2][1] += E2[1];
        electricField[pair.dipole2][2] += E2[2];
    }
}

double DrudeSCFOpenMM::updateDrudePositions(model::MCState& state,
                                            const std::vector<DrudeParticle>& particles,
                                            const std::vector<Vec3>& electricField,
                                            double damping,
                                            double maxStep,
                                            double hardWall) const {
    double maxDisplacement = 0.0;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];
        
        // Freeze particles with zero polarizability
        if (particle.polarizability < 1e-14) {
            drude.x = parent.x;
            drude.y = parent.y;
            drude.z = parent.z;
            continue;
        }
        
        // Target displacement: d = q*E/k
        double targetX = particle.charge * electricField[i][0] / particle.kSpring;
        double targetY = particle.charge * electricField[i][1] / particle.kSpring;
        double targetZ = particle.charge * electricField[i][2] / particle.kSpring;
        
        // Current displacement
        double oldX = drude.x - parent.x;
        double oldY = drude.y - parent.y;
        double oldZ = drude.z - parent.z;
        
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(oldX, oldY, oldZ, box);
        
        // Mixed update with damping
        double newX = (1.0 - damping) * oldX + damping * targetX;
        double newY = (1.0 - damping) * oldY + damping * targetY;
        double newZ = (1.0 - damping) * oldZ + damping * targetZ;
        
        // Apply step size limit
        double stepX = newX - oldX;
        double stepY = newY - oldY;
        double stepZ = newZ - oldZ;
        double step2 = stepX*stepX + stepY*stepY + stepZ*stepZ;
        
        if (step2 > maxStep * maxStep) {
            double scale = maxStep / std::sqrt(step2);
            newX = oldX + stepX * scale;
            newY = oldY + stepY * scale;
            newZ = oldZ + stepZ * scale;
        }
        
        // Apply hard wall constraint
        double r2 = newX*newX + newY*newY + newZ*newZ;
        if (hardWall < 1e8 && r2 > hardWall * hardWall) {
            double scale = hardWall / std::sqrt(r2);
            newX *= scale;
            newY *= scale;
            newZ *= scale;
        }
        
        // Update position
        drude.x = parent.x + newX;
        drude.y = parent.y + newY;
        drude.z = parent.z + newZ;
        
        // Track maximum displacement
        double change2 = (newX - oldX) * (newX - oldX) +
                        (newY - oldY) * (newY - oldY) +
                        (newZ - oldZ) * (newZ - oldZ);
        maxDisplacement = std::max(maxDisplacement, std::sqrt(change2));
    }
    
    return maxDisplacement;
}

double DrudeSCFOpenMM::calculateSpringEnergy(const model::MCState& state,
                                             const std::vector<DrudeParticle>& particles) const {
    double energy = 0.0;
    
    for (const auto& p : particles) {
        const auto& drude = state.atoms[p.drudeIndex];
        const auto& parent = state.atoms[p.parentIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        energy += 0.5 * p.kSpring * r2;
    }
    
    return energy;
}

bool DrudeSCFOpenMM::inSameMolecule(int atom1, int atom2, const model::MCState& state) const {
    // Simple implementation - check if atoms belong to same residue
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (i >= static_cast<int>(state.residues.size())) {
            break;
        }
        
        const auto& res = state.residues[i];
        int start = res.atomStart;
        int end = start + res.atomCount;
        
        bool atom1InRes = (atom1 >= start && atom1 < end);
        bool atom2InRes = (atom2 >= start && atom2 < end);
        
        if (atom1InRes && atom2InRes) {
            return true;
        }
    }
    
    return false;
}

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc