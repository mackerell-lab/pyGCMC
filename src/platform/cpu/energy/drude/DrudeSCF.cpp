/**
 * @file DrudeSCF.cpp
 * @brief Self-Consistent Field (SCF) optimizer implementation
 */

#include "DrudeSCF.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

void DrudeSCF::buildActiveAtomMask(const model::MCState& state, std::vector<char>& mask) {
    const int nAtoms = std::max(0, state.activeAtomCount);
    mask.assign(static_cast<size_t>(nAtoms), 1);
    if (state.residues.empty()) {
        return;
    }

    // Default: all atoms [0..activeAtomCount) are active.
    // Mask out atoms belonging to inactive residues (GCMC deletion keeps atom arrays un-compacted),
    // while leaving atoms not covered by any residue range active (tests may model external charges this way).
    for (const auto& res : state.residues) {
        if (res.active) {
            continue;
        }
        if (res.atomCount <= 0) {
            continue;
        }
        const int start = std::max(0, res.atomStart);
        const int end = std::min(nAtoms, res.atomStart + res.atomCount);
        for (int j = start; j < end; ++j) {
            mask[static_cast<size_t>(j)] = 0;
        }
    }
}

void DrudeSCF::buildActiveAtomIndices(const std::vector<char>& mask, std::vector<int>& indices) {
    indices.clear();
    indices.reserve(mask.size());
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
            indices.push_back(static_cast<int>(i));
        }
    }
}

bool DrudeSCF::isActiveAtom(int atomIndex, const std::vector<char>& mask) {
    if (atomIndex < 0) {
        return false;
    }
    const size_t idx = static_cast<size_t>(atomIndex);
    if (idx >= mask.size()) {
        return false;
    }
    return mask[idx] != 0;
}

bool DrudeSCF::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    if (particles.empty()) return true;

    std::vector<char> activeAtomMask;
    std::vector<int> activeAtoms;
    buildActiveAtomMask(state, activeAtomMask);
    buildActiveAtomIndices(activeAtomMask, activeAtoms);
    
    std::vector<Vec3> electricField(particles.size(), {0.0, 0.0, 0.0});
    std::vector<Vec3> drudeForces(particles.size(), {0.0, 0.0, 0.0});
    
    // Initialize Drude positions near parents if needed
    for (const auto& particle : particles) {
        if (!isActiveAtom(particle.parentIndex, activeAtomMask) ||
            !isActiveAtom(particle.drudeIndex, activeAtomMask)) {
            continue;
        }
        
        const auto& parent = state.atoms[particle.parentIndex];
        auto& drude = state.atoms[particle.drudeIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double dist2 = dx*dx + dy*dy + dz*dz;
        if (params.enableHardWall && dist2 > params.maxDrudeDistance * params.maxDrudeDistance) {
            // Reset Drude to parent position
            drude.x = parent.x;
            drude.y = parent.y;
            drude.z = parent.z;
        }
    }
    
    // SCF iteration
    double dampingFactor = params.dampingFactor;
    m_lastIterationCount = 0;
    
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        m_lastIterationCount = iter + 1;
        // Calculate electric field at each Drude particle
        std::fill(electricField.begin(), electricField.end(), Vec3{0.0, 0.0, 0.0});
        calculateElectricField(state, particles, screenedPairs, electricField, activeAtoms, activeAtomMask);
        
        // Calculate forces on Drude particles
        calculateDrudeForces(state, particles, electricField, drudeForces, activeAtomMask);
        
        // Check convergence with hard wall consideration
        double maxForce = 0.0;
        bool anyAtHardWall = false;
        
        for (size_t i = 0; i < particles.size(); ++i) {
            const auto& particle = particles[i];
            if (!isActiveAtom(particle.parentIndex, activeAtomMask) ||
                !isActiveAtom(particle.drudeIndex, activeAtomMask)) {
                continue;
            }
            const auto& drude = state.atoms[particle.drudeIndex];
            const auto& parent = state.atoms[particle.parentIndex];
            
            // Check if at hard wall
            double dx = drude.x - parent.x;
            double dy = drude.y - parent.y;
            double dz = drude.z - parent.z;
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(dx, dy, dz, box);
            double dist2 = dx*dx + dy*dy + dz*dz;
            
            if (params.enableHardWall && dist2 > 0.98 * params.maxDrudeDistance * params.maxDrudeDistance) {
                anyAtHardWall = true;
            }
            
            double f2 = drudeForces[i][0]*drudeForces[i][0] + 
                       drudeForces[i][1]*drudeForces[i][1] + 
                       drudeForces[i][2]*drudeForces[i][2];
            maxForce = std::max(maxForce, f2);
        }
        maxForce = std::sqrt(maxForce);
        
        // Relax convergence criteria if at hard wall
        double effectiveTolerance = anyAtHardWall ? params.tolerance * 10.0 : params.tolerance;
        
        if (maxForce < effectiveTolerance) {
            return true;  // Converged
        }
        
        // Debug output - disabled for now
        // if (iter < 5 || iter % 10 == 0) {
        //     std::cerr << "SCF iter " << iter << ": maxForce = " << maxForce 
        //               << ", tolerance = " << params.tolerance << std::endl;
        // }
        
        // Update Drude positions
        updateDrudePositions(
            state, particles, electricField, dampingFactor, 
            params.enableHardWall ? params.maxDrudeDistance : 1e10,  // Large value when hard wall disabled
            activeAtomMask
        );
        
        // Adaptive damping - disabled for now to debug
        // if (iter > 5 && maxDisplacement > 0.001) {
        //     dampingFactor = std::max(0.1, dampingFactor * 0.95);
        // }
    }
    
    return false;  // Did not converge
}

void DrudeSCF::calculateElectricField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    std::vector<Vec3>& electricField,
    const std::vector<int>& activeAtoms,
    const std::vector<char>& activeAtomMask
) const {
    // Field from external charges
    calculateExternalField(state, particles, electricField, activeAtoms);
    
    // Field from induced dipoles (other Drude particles)
    calculateInducedField(state, particles, screenedPairs, electricField, activeAtomMask);
}

void DrudeSCF::calculateExternalField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& electricField,
    const std::vector<int>& activeAtoms
) const {
    // Loop over all Drude particles
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        
        // Bounds check for drude index
        if (particle.drudeIndex < 0 || particle.drudeIndex >= state.activeAtomCount ||
            particle.parentIndex < 0 || particle.parentIndex >= state.activeAtomCount) {
            continue;
        }
        
        const auto& drudeAtom = state.atoms[particle.drudeIndex];
        
        Vec3 field = {0.0, 0.0, 0.0};
        
        // Loop over active atoms (skip inactive/ghost atoms).
        for (int j : activeAtoms) {
            // Skip the Drude particle itself
            if (j == particle.drudeIndex) continue;
            
            // Skip the parent atom (Drude-parent interaction is handled by spring)
            if (j == particle.parentIndex) continue;
            
            // Skip intramolecular interactions
            // Drude particles should NOT interact with other atoms in the same molecule
            if (inSameMolecule(particle.drudeIndex, j, state)) continue;
            
            const auto& atom = state.atoms[j];
            
            // Calculate distance
            double dx = drudeAtom.x - atom.x;  // Fixed: vector from source to field point
            double dy = drudeAtom.y - atom.y;
            double dz = drudeAtom.z - atom.z;
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;  // Skip if too close
            
            double r = std::sqrt(r2);
            double r3 = r2 * r;
            
            // Electric field: E = k * q / r^2 * r_hat
            // r_hat points from source to field point
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

void DrudeSCF::calculateInducedField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    std::vector<Vec3>& electricField,
    const std::vector<char>& activeAtomMask
) const {
    // External field already includes *unscreened* Coulomb contributions from all atoms
    // (including other Drude and parent atoms). Here we apply a *correction* for
    // the configured screened dipole pairs: add (screened - unscreened) field.
    //
    // For a screened Coulomb energy U = k*q1*q2*S1(u)/r with u = a*r/α_eff,
    // the corresponding force/field damping factor is:
    //   D(u) = S1(u) - u*S1'(u) = 1 - (1 + u + u^2/2) * exp(-u)
    // This fixes the historical missing (u^2/2) term and avoids double counting.
    const std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};

    for (const auto& pair : screenedPairs) {
        if (pair.dipole1 < 0 || static_cast<size_t>(pair.dipole1) >= particles.size() ||
            pair.dipole2 < 0 || static_cast<size_t>(pair.dipole2) >= particles.size()) {
            continue;
        }

        const auto& dipole1 = particles[pair.dipole1];
        const auto& dipole2 = particles[pair.dipole2];

        if (!isActiveAtom(dipole1.drudeIndex, activeAtomMask) ||
            !isActiveAtom(dipole1.parentIndex, activeAtomMask) ||
            !isActiveAtom(dipole2.drudeIndex, activeAtomMask) ||
            !isActiveAtom(dipole2.parentIndex, activeAtomMask)) {
            continue;
        }

        const double alphaEff = std::pow(dipole1.polarizability * dipole2.polarizability, 1.0 / 6.0);
        if (alphaEff < 1e-14) {
            continue;
        }

        auto applyCorrectionFromSource = [&](int fieldDipoleIndex, int fieldAtomIndex, int sourceAtomIndex) {
            // Match calculateExternalField() skip semantics: if the base unscreened term
            // was skipped, do not attempt to "correct" it.
            if (!isActiveAtom(fieldAtomIndex, activeAtomMask) ||
                !isActiveAtom(sourceAtomIndex, activeAtomMask)) {
                return;
            }
            if (inSameMolecule(fieldAtomIndex, sourceAtomIndex, state)) {
                return;
            }

            const auto& fieldAtom = state.atoms[fieldAtomIndex];
            const auto& sourceAtom = state.atoms[sourceAtomIndex];

            if (std::abs(sourceAtom.charge) < 1e-14) {
                return;
            }

            double dx = fieldAtom.x - sourceAtom.x;
            double dy = fieldAtom.y - sourceAtom.y;
            double dz = fieldAtom.z - sourceAtom.z;
            applyPBC(dx, dy, dz, box);

            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < 1e-12) {
                return;
            }

            const double r = std::sqrt(r2);
            const double invR3 = 1.0 / (r2 * r);

            // Unscreened field contribution from the source charge.
            const double unscreenedFactor = DrudeConstants::ONE_4PI_EPS0 * sourceAtom.charge * invR3;

            // Compute D(u) for screened Coulomb field (screened - unscreened correction).
            double damping = 1.0;
            if (pair.thole != 0.0) {
                const double u = pair.thole * r / alphaEff;
                if (u <= 50.0) {
                    const double expu = std::exp(-u);
                    damping = 1.0 - expu * (1.0 + u + 0.5 * u * u);
                }
            }

            const double deltaFactor = (damping - 1.0) * unscreenedFactor;
            electricField[fieldDipoleIndex][0] += deltaFactor * dx;
            electricField[fieldDipoleIndex][1] += deltaFactor * dy;
            electricField[fieldDipoleIndex][2] += deltaFactor * dz;
        };

        // Field at Drude of dipole1: correct interactions with parent2 and drude2.
        applyCorrectionFromSource(pair.dipole1, dipole1.drudeIndex, dipole2.parentIndex);
        applyCorrectionFromSource(pair.dipole1, dipole1.drudeIndex, dipole2.drudeIndex);

        // Field at Drude of dipole2: correct interactions with parent1 and drude1.
        applyCorrectionFromSource(pair.dipole2, dipole2.drudeIndex, dipole1.parentIndex);
        applyCorrectionFromSource(pair.dipole2, dipole2.drudeIndex, dipole1.drudeIndex);
    }
}

double DrudeSCF::updateDrudePositions(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& electricField,
    double dampingFactor,
    double maxDrudeDistance,
    const std::vector<char>& activeAtomMask
) const {
    double maxDisplacement = 0.0;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];

        if (!isActiveAtom(particle.parentIndex, activeAtomMask) ||
            !isActiveAtom(particle.drudeIndex, activeAtomMask)) {
            continue;
        }
        
        auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];
        
        // New equilibrium position: r_drude = r_parent + q * E / k
        // For negative charge Drude, this gives displacement opposite to E field
        double dx_new = particle.charge * electricField[i][0] / particle.kSpring;
        double dy_new = particle.charge * electricField[i][1] / particle.kSpring;
        double dz_new = particle.charge * electricField[i][2] / particle.kSpring;
        
        // Don't limit dx_new here - let SCF find the true equilibrium
        // Hard wall will be applied after convergence if needed
        
        // Current displacement
        double dx_old = drude.x - parent.x;
        double dy_old = drude.y - parent.y;
        double dz_old = drude.z - parent.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx_old, dy_old, dz_old, box);
        
        // Apply damping
        double dx_update = (1.0 - dampingFactor) * dx_old + dampingFactor * dx_new;
        double dy_update = (1.0 - dampingFactor) * dy_old + dampingFactor * dy_new;
        double dz_update = (1.0 - dampingFactor) * dz_old + dampingFactor * dz_new;
        
        // Update position with hard wall constraint (if enabled)
        // Apply hard wall only to the final position, not the target
        double disp2_update = dx_update*dx_update + dy_update*dy_update + dz_update*dz_update;
        if (maxDrudeDistance < 1.0 && disp2_update > maxDrudeDistance * maxDrudeDistance) {
            double scale = maxDrudeDistance / std::sqrt(disp2_update);
            dx_update *= scale;
            dy_update *= scale;
            dz_update *= scale;
        }
        
        // Update position
        drude.x = parent.x + dx_update;
        drude.y = parent.y + dy_update;
        drude.z = parent.z + dz_update;
        
        // Track maximum displacement
        double change2 = (dx_update - dx_old) * (dx_update - dx_old) +
                        (dy_update - dy_old) * (dy_update - dy_old) +
                        (dz_update - dz_old) * (dz_update - dz_old);
        maxDisplacement = std::max(maxDisplacement, std::sqrt(change2));
    }
    
    return maxDisplacement;
}

void DrudeSCF::calculateDrudeForces(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& electricField,
    std::vector<Vec3>& forces,
    const std::vector<char>& activeAtomMask
) const {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        if (!isActiveAtom(particle.parentIndex, activeAtomMask) ||
            !isActiveAtom(particle.drudeIndex, activeAtomMask)) {
            forces[i] = {0.0, 0.0, 0.0};
            continue;
        }
        const auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];
        
        // Spring force: F_spring = -k * (r_drude - r_parent)
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        Vec3 springForce = {-particle.kSpring * dx,
                           -particle.kSpring * dy,
                           -particle.kSpring * dz};
        
        // Electric force: F_elec = q * E
        Vec3 electricForce = {particle.charge * electricField[i][0],
                             particle.charge * electricField[i][1],
                             particle.charge * electricField[i][2]};
        
        // Total force
        forces[i][0] = springForce[0] + electricForce[0];
        forces[i][1] = springForce[1] + electricForce[1];
        forces[i][2] = springForce[2] + electricForce[2];
    }
}

void DrudeSCF::applyPBC(double& dx, double& dy, double& dz, const std::array<double, 3>& box) const {
    if (box[0] > 0) dx -= box[0] * std::round(dx / box[0]);
    if (box[1] > 0) dy -= box[1] * std::round(dy / box[1]);
    if (box[2] > 0) dz -= box[2] * std::round(dz / box[2]);
}

bool DrudeSCF::inSameMolecule(int atom1, int atom2, const model::MCState& state) const {
    // Check if we have residues defined
    if (state.activeResidueCount <= 0 || state.residues.empty()) {
        return false;  // No residues defined, assume all atoms are in different molecules
    }
    
    // Find which residue each atom belongs to
    for (int i = 0; i < state.activeResidueCount; ++i) {
        // Bounds check
        if (i >= static_cast<int>(state.residues.size())) {
            break;
        }
        
        const auto& res = state.residues[i];
        if (!res.active) {
            continue;
        }
        int start = res.atomStart;
        int end = start + res.atomCount;
        
        bool atom1InRes = (atom1 >= start && atom1 < end);
        bool atom2InRes = (atom2 >= start && atom2 < end);
        
        if (atom1InRes && atom2InRes) {
            return true;  // Both atoms in same residue
        }
    }
    
    return false;  // Atoms in different residues
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
