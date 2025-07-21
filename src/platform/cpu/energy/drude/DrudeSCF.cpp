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

bool DrudeSCF::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    if (particles.empty()) return true;
    
    std::vector<Vec3> electricField(particles.size(), {0.0, 0.0, 0.0});
    std::vector<Vec3> drudeForces(particles.size(), {0.0, 0.0, 0.0});
    
    // Initialize Drude positions near parents if needed
    for (const auto& particle : particles) {
        const auto& parent = state.atoms[particle.parentIndex];
        auto& drude = state.atoms[particle.drudeIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double dist2 = dx*dx + dy*dy + dz*dz;
        if (dist2 > params.maxDrudeDistance * params.maxDrudeDistance) {
            // Reset Drude to parent position
            drude.x = parent.x;
            drude.y = parent.y;
            drude.z = parent.z;
        }
    }
    
    // SCF iteration
    double dampingFactor = params.dampingFactor;
    
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        // Calculate electric field at each Drude particle
        std::fill(electricField.begin(), electricField.end(), Vec3{0.0, 0.0, 0.0});
        calculateElectricField(state, particles, screenedPairs, electricField);
        
        // Calculate forces on Drude particles
        calculateDrudeForces(state, particles, electricField, drudeForces);
        
        // Check convergence
        double maxForce = 0.0;
        for (const auto& force : drudeForces) {
            double f2 = force[0]*force[0] + force[1]*force[1] + force[2]*force[2];
            maxForce = std::max(maxForce, f2);
        }
        maxForce = std::sqrt(maxForce);
        
        if (maxForce < params.tolerance) {
            return true;  // Converged
        }
        
        // Update Drude positions
        double maxDisplacement = updateDrudePositions(
            state, particles, electricField, dampingFactor, params.maxDrudeDistance
        );
        
        // Adaptive damping
        if (iter > 5 && maxDisplacement > 0.001) {
            dampingFactor = std::max(0.1, dampingFactor * 0.95);
        }
    }
    
    return false;  // Did not converge
}

void DrudeSCF::calculateElectricField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    std::vector<Vec3>& electricField
) const {
    // Field from external charges
    calculateExternalField(state, particles, electricField);
    
    // Field from induced dipoles (other Drude particles)
    calculateInducedField(state, particles, screenedPairs, electricField);
}

void DrudeSCF::calculateExternalField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& electricField
) const {
    // Loop over all Drude particles
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& drudeAtom = state.atoms[particle.drudeIndex];
        
        Vec3 field = {0.0, 0.0, 0.0};
        
        // Loop over all atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip the Drude particle itself
            if (j == particle.drudeIndex) continue;
            
            // Skip if in the same molecule (intramolecular exclusion)
            if (inSameMolecule(particle.drudeIndex, j, state)) continue;
            
            const auto& atom = state.atoms[j];
            
            // Calculate distance
            double dx = atom.x - drudeAtom.x;
            double dy = atom.y - drudeAtom.y;
            double dz = atom.z - drudeAtom.z;
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

void DrudeSCF::calculateInducedField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    std::vector<Vec3>& electricField
) const {
    // Field from Drude-Drude interactions with Thole screening
    for (const auto& pair : screenedPairs) {
        const auto& particle1 = particles[pair.dipole1];
        const auto& particle2 = particles[pair.dipole2];
        
        // Bounds check
        if (particle1.drudeIndex < 0 || particle1.drudeIndex >= state.activeAtomCount ||
            particle2.drudeIndex < 0 || particle2.drudeIndex >= state.activeAtomCount) {
            continue;
        }
        
        const auto& drude1 = state.atoms[particle1.drudeIndex];
        const auto& drude2 = state.atoms[particle2.drudeIndex];
        
        // Calculate distance
        double dx = drude2.x - drude1.x;
        double dy = drude2.y - drude1.y;
        double dz = drude2.z - drude1.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-12) continue;
        
        double r = std::sqrt(r2);
        
        // Calculate Thole screening and its derivative
        double alpha_ij = std::pow(particle1.polarizability * particle2.polarizability, 1.0/6.0);
        double u = r / alpha_ij;
        double thu = pair.thole * u;
        double exp_thu = std::exp(-thu);
        
        // Thole damping function for dipole field
        double damping = 1.0 - exp_thu * (1.0 + thu);
        
        // Electric field from dipole 2 at dipole 1
        double factor = DrudeConstants::ONE_4PI_EPS0 * particle2.charge * damping / (r2 * r);
        
        electricField[pair.dipole1][0] += factor * dx;
        electricField[pair.dipole1][1] += factor * dy;
        electricField[pair.dipole1][2] += factor * dz;
        
        // Electric field from dipole 1 at dipole 2 (Newton's third law)
        factor = DrudeConstants::ONE_4PI_EPS0 * particle1.charge * damping / (r2 * r);
        
        electricField[pair.dipole2][0] -= factor * dx;
        electricField[pair.dipole2][1] -= factor * dy;
        electricField[pair.dipole2][2] -= factor * dz;
    }
}

double DrudeSCF::updateDrudePositions(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& electricField,
    double dampingFactor,
    double maxDrudeDistance
) const {
    double maxDisplacement = 0.0;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];
        
        // New equilibrium position: r_drude = r_parent + q * E / k
        double dx_new = particle.charge * electricField[i][0] / particle.kSpring;
        double dy_new = particle.charge * electricField[i][1] / particle.kSpring;
        double dz_new = particle.charge * electricField[i][2] / particle.kSpring;
        
        // Limit displacement magnitude
        double disp2 = dx_new*dx_new + dy_new*dy_new + dz_new*dz_new;
        if (disp2 > maxDrudeDistance * maxDrudeDistance) {
            double scale = maxDrudeDistance / std::sqrt(disp2);
            dx_new *= scale;
            dy_new *= scale;
            dz_new *= scale;
        }
        
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
    std::vector<Vec3>& forces
) const {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
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
    dx -= box[0] * std::round(dx / box[0]);
    dy -= box[1] * std::round(dy / box[1]);
    dz -= box[2] * std::round(dz / box[2]);
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