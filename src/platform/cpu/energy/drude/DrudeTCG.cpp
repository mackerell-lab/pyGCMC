/**
 * @file DrudeTCG.cpp
 * @brief Truncated Conjugate Gradient implementation
 */

#include "DrudeTCG.hpp"
#include "../common/EnergyConstants.hpp"
#include <cmath>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

bool DrudeTCG::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    // Mark screenedPairs as unused for now
    (void)screenedPairs;
    (void)params;
    
    const size_t nParticles = particles.size();
    if (nParticles == 0) return true;
    
    // Allocate working arrays if needed
    if (residuals_.size() != nParticles) {
        residuals_.resize(nParticles);
        directions_.resize(nParticles);
        Ap_.resize(nParticles);
        oldResiduals_.resize(nParticles);
    }
    
    // For Drude model, we optimize positions directly
    // The force balance equation is: k*d = q*E
    // We solve this using CG on the positions
    
    // Step 1: Compute permanent fields E0 at parent positions
    std::vector<Vec3> E0(nParticles);
    computeFieldsAtParents(state, particles, E0);
    
    // Step 2: Initialize CG for position optimization
    // We solve: (k*I - q*q*T) * d = q*E0
    // where d is the Drude displacement from parent
    
    // Initial guess: d = 0 (Drude at parent position)
    std::vector<Vec3> displacements(nParticles, {0.0, 0.0, 0.0});
    
    // Initial residual: r = q*E0/k (target displacement)
    // Initial search direction: p = r
    for (size_t i = 0; i < nParticles; ++i) {
        const auto& particle = particles[i];
        double factor = particle.charge / particle.kSpring;
        residuals_[i][0] = factor * E0[i][0];
        residuals_[i][1] = factor * E0[i][1];
        residuals_[i][2] = factor * E0[i][2];
        directions_[i] = residuals_[i];
    }
    
    // Step 3: CG iterations for position optimization
    double rsOld = 0.0;
    for (size_t i = 0; i < nParticles; ++i) {
        rsOld += residuals_[i][0] * residuals_[i][0] +
                 residuals_[i][1] * residuals_[i][1] +
                 residuals_[i][2] * residuals_[i][2];
    }
    
    for (int iter = 0; iter < tcgIterations_; ++iter) {
        // Update Drude positions for field calculation
        updateDrudePositions(state, particles, displacements);
        
        // Compute A*p where A = I - (q/k)*T
        // First term: p
        for (size_t i = 0; i < nParticles; ++i) {
            Ap_[i] = directions_[i];
        }
        
        // Second term: -(q/k)*T*p
        // T*p means fields at Drude positions due to displacements p
        std::vector<Vec3> inducedFields(nParticles);
        computeInducedFieldsFromDisplacements(state, particles, directions_, inducedFields);
        
        for (size_t i = 0; i < nParticles; ++i) {
            const auto& particle = particles[i];
            double factor = particle.charge / particle.kSpring;
            Ap_[i][0] -= factor * inducedFields[i][0];
            Ap_[i][1] -= factor * inducedFields[i][1];
            Ap_[i][2] -= factor * inducedFields[i][2];
        }
        
        // Compute step size: alpha = rsOld / (p^T * Ap)
        double pAp = 0.0;
        for (size_t i = 0; i < nParticles; ++i) {
            pAp += directions_[i][0] * Ap_[i][0] +
                   directions_[i][1] * Ap_[i][1] +
                   directions_[i][2] * Ap_[i][2];
        }
        
        if (std::abs(pAp) < 1e-10) break;
        
        double alpha = rsOld / pAp;
        
        // Update displacements and residuals
        for (size_t i = 0; i < nParticles; ++i) {
            // d = d + alpha * p
            displacements[i][0] += alpha * directions_[i][0];
            displacements[i][1] += alpha * directions_[i][1];
            displacements[i][2] += alpha * directions_[i][2];
            
            // r = r - alpha * Ap
            residuals_[i][0] -= alpha * Ap_[i][0];
            residuals_[i][1] -= alpha * Ap_[i][1];
            residuals_[i][2] -= alpha * Ap_[i][2];
        }
        
        // Compute new rsNew
        double rsNew = 0.0;
        for (size_t i = 0; i < nParticles; ++i) {
            rsNew += residuals_[i][0] * residuals_[i][0] +
                     residuals_[i][1] * residuals_[i][1] +
                     residuals_[i][2] * residuals_[i][2];
        }
        
        // Update search direction: p = r + beta * p
        double beta = rsNew / rsOld;
        for (size_t i = 0; i < nParticles; ++i) {
            directions_[i][0] = residuals_[i][0] + beta * directions_[i][0];
            directions_[i][1] = residuals_[i][1] + beta * directions_[i][1];
            directions_[i][2] = residuals_[i][2] + beta * directions_[i][2];
        }
        
        rsOld = rsNew;
    }
    
    // Step 4: Apply final displacements to update Drude positions
    updateDrudePositions(state, particles, displacements);
    
    return true;
}

void DrudeTCG::computeFieldsAtParents(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& fields
) {
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};
    
    // Initialize fields to zero
    for (auto& field : fields) {
        field[0] = field[1] = field[2] = 0.0;
    }
    
    // Compute field at each parent position from all charges
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parent = state.atoms[particle.parentIndex];
        
        // Sum over all atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self (parent)
            if (j == particle.parentIndex) continue;
            
            // Skip Drude (will be at parent initially)
            if (j == particle.drudeIndex) continue;
            
            const auto& atom = state.atoms[j];
            
            // Skip if no charge
            if (std::abs(atom.charge) < 1e-6) continue;
            
            // Compute distance with PBC
            double dx = parent.x - atom.x;
            double dy = parent.y - atom.y;
            double dz = parent.z - atom.z;
            
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
            
            fields[i][0] += fieldMag * dx;
            fields[i][1] += fieldMag * dy;
            fields[i][2] += fieldMag * dz;
        }
    }
}

void DrudeTCG::updateDrudePositions(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& displacements
) {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parent = state.atoms[particle.parentIndex];
        auto& drude = state.atoms[particle.drudeIndex];
        
        // Update Drude position
        drude.x = parent.x + displacements[i][0];
        drude.y = parent.y + displacements[i][1];
        drude.z = parent.z + displacements[i][2];
    }
}

void DrudeTCG::computeInducedFieldsFromDisplacements(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& displacements,
    std::vector<Vec3>& fields
) {
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};
    
    // Initialize
    for (auto& field : fields) {
        field[0] = field[1] = field[2] = 0.0;
    }
    
    // Compute fields at Drude positions due to other Drude charges
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle_i = particles[i];
        const auto& parent_i = state.atoms[particle_i.parentIndex];
        
        // Position of Drude i
        double xi = parent_i.x + displacements[i][0];
        double yi = parent_i.y + displacements[i][1];
        double zi = parent_i.z + displacements[i][2];
        
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;
            
            const auto& particle_j = particles[j];
            const auto& parent_j = state.atoms[particle_j.parentIndex];
            
            // Position of Drude j
            double xj = parent_j.x + displacements[j][0];
            double yj = parent_j.y + displacements[j][1];
            double zj = parent_j.z + displacements[j][2];
            
            // Distance with PBC
            double dx = xi - xj;
            double dy = yi - yj;
            double dz = zi - zj;
            
            if (dx > halfBox[0]) dx -= box[0];
            if (dx < -halfBox[0]) dx += box[0];
            if (dy > halfBox[1]) dy -= box[1];
            if (dy < -halfBox[1]) dy += box[1];
            if (dz > halfBox[2]) dz -= box[2];
            if (dz < -halfBox[2]) dz += box[2];
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > cutoff2 || r2 < 1e-10) continue;
            
            // E = k*q/r^2 * r_hat
            double r = std::sqrt(r2);
            double fieldMag = DrudeConstants::ONE_4PI_EPS0 * particle_j.charge / (r2 * r);
            
            fields[i][0] += fieldMag * dx;
            fields[i][1] += fieldMag * dy;
            fields[i][2] += fieldMag * dz;
        }
    }
}

void DrudeTCG::computeFields(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& fields
) {
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};
    
    // Initialize fields to zero
    for (auto& field : fields) {
        field[0] = field[1] = field[2] = 0.0;
    }
    
    // Compute field at each Drude position from all charges
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& drude = state.atoms[particle.drudeIndex];
        
        // Sum over all atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self
            if (j == particle.drudeIndex) continue;
            
            const auto& atom = state.atoms[j];
            
            // Skip if no charge
            if (std::abs(atom.charge) < 1e-6) continue;
            
            // Skip parent-Drude interaction
            if (j == particle.parentIndex) continue;
            
            // Compute distance with PBC
            double dx = drude.x - atom.x;
            double dy = drude.y - atom.y;
            double dz = drude.z - atom.z;
            
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
            
            fields[i][0] += fieldMag * dx;
            fields[i][1] += fieldMag * dy;
            fields[i][2] += fieldMag * dz;
        }
    }
}

bool DrudeTCG::isExcluded(int atom1, int atom2, const model::MCState& state) {
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
    
    // Excluded if in same residue
    return (res1 >= 0 && res1 == res2);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc