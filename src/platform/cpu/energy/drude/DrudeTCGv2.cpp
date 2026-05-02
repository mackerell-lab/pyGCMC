/**
 * @file DrudeTCGv2.cpp
 * @brief Improved TCG implementation with better accuracy
 */

#include "DrudeTCGv2.hpp"
#include "../common/EnergyConstants.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

bool DrudeTCGv2::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    (void)screenedPairs;
    (void)params;

    const size_t nParticles = particles.size();
    if (nParticles == 0) return true;

    // Enable debug for small systems
    if (nParticles <= 2) {
        debugMode_ = true;
    }

    // Allocate working arrays
    if (residuals_.size() != nParticles) {
        residuals_.resize(nParticles);
        directions_.resize(nParticles);
        Ap_.resize(nParticles);
        z_.resize(nParticles);
        precond_.resize(nParticles);
    }

    // Build diagonal preconditioner
    if (usePreconditioning_) {
        buildPreconditioner(particles);
    }

    // Step 1: Compute permanent fields E0 at parent positions
    std::vector<Vec3> E0(nParticles);
    computeFieldsAtParents(state, particles, E0);

    // Step 2: Initialize for CG iteration
    // We solve: A*d = b where A = I - (q/k)*T and b = (q/k)*E0

    // Initial guess: d = 0 (Drude at parent position)
    std::vector<Vec3> displacements(nParticles, {0.0, 0.0, 0.0});

    // Initial residual: r = b = (q/k)*E0
    for (size_t i = 0; i < nParticles; ++i) {
        const auto& particle = particles[i];
        double factor = particle.charge / particle.kSpring;
        residuals_[i][0] = factor * E0[i][0];
        residuals_[i][1] = factor * E0[i][1];
        residuals_[i][2] = factor * E0[i][2];
    }

    // Apply preconditioner: z = M^{-1} * r
    if (usePreconditioning_) {
        applyPreconditioner(residuals_, z_);
        // Initial direction: p = z
        directions_ = z_;
    } else {
        // Initial direction: p = r
        directions_ = residuals_;
    }

    // Compute initial residual norm for preconditioned CG
    double rzOld = 0.0;
    if (usePreconditioning_) {
        for (size_t i = 0; i < nParticles; ++i) {
            rzOld += residuals_[i][0] * z_[i][0] +
                     residuals_[i][1] * z_[i][1] +
                     residuals_[i][2] * z_[i][2];
        }
    } else {
        for (size_t i = 0; i < nParticles; ++i) {
            rzOld += residuals_[i][0] * residuals_[i][0] +
                     residuals_[i][1] * residuals_[i][1] +
                     residuals_[i][2] * residuals_[i][2];
        }
    }

    if (debugMode_) {
        std::cout << "\nTCGv2 Debug Output:" << std::endl;
        std::cout << "Number of particles: " << nParticles << std::endl;
        std::cout << "Initial fields E0:" << std::endl;
        for (size_t i = 0; i < nParticles; ++i) {
            std::cout << "  Particle " << i << ": ("
                      << E0[i][0] << ", " << E0[i][1] << ", " << E0[i][2]
                      << ") mag=" << std::sqrt(E0[i][0]*E0[i][0] + E0[i][1]*E0[i][1] + E0[i][2]*E0[i][2])
                      << std::endl;
        }
        std::cout << "Initial residual norm = " << std::sqrt(rzOld) << std::endl;
    }

    // Step 3: CG iterations
    for (int iter = 0; iter < tcgIterations_; ++iter) {
        // Update Drude positions for correct field calculation
        updateDrudePositions(state, particles, displacements);

        // Compute A*p
        applyMatrix(state, particles, directions_, Ap_);

        // Compute step size: alpha = (r·z) / (p·Ap)
        double pAp = 0.0;
        for (size_t i = 0; i < nParticles; ++i) {
            pAp += directions_[i][0] * Ap_[i][0] +
                   directions_[i][1] * Ap_[i][1] +
                   directions_[i][2] * Ap_[i][2];
        }

        if (std::abs(pAp) < 1e-10) {
            if (debugMode_) {
                std::cout << "TCGv2: Converged at iteration " << iter << std::endl;
            }
            break;
        }

        double alpha = rzOld / pAp;

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

        // Apply preconditioner to new residual
        if (usePreconditioning_) {
            applyPreconditioner(residuals_, z_);
        }

        // Compute new (r·z)
        double rzNew = 0.0;
        if (usePreconditioning_) {
            for (size_t i = 0; i < nParticles; ++i) {
                rzNew += residuals_[i][0] * z_[i][0] +
                         residuals_[i][1] * z_[i][1] +
                         residuals_[i][2] * z_[i][2];
            }
        } else {
            for (size_t i = 0; i < nParticles; ++i) {
                rzNew += residuals_[i][0] * residuals_[i][0] +
                         residuals_[i][1] * residuals_[i][1] +
                         residuals_[i][2] * residuals_[i][2];
            }
        }

        if (debugMode_) {
            std::cout << "TCGv2: Iteration " << iter + 1
                      << ", residual norm = " << std::sqrt(rzNew) << std::endl;
        }

        // Update search direction: p = z + beta * p (or p = r + beta * p)
        double beta = rzNew / rzOld;
        if (usePreconditioning_) {
            for (size_t i = 0; i < nParticles; ++i) {
                directions_[i][0] = z_[i][0] + beta * directions_[i][0];
                directions_[i][1] = z_[i][1] + beta * directions_[i][1];
                directions_[i][2] = z_[i][2] + beta * directions_[i][2];
            }
        } else {
            for (size_t i = 0; i < nParticles; ++i) {
                directions_[i][0] = residuals_[i][0] + beta * directions_[i][0];
                directions_[i][1] = residuals_[i][1] + beta * directions_[i][1];
                directions_[i][2] = residuals_[i][2] + beta * directions_[i][2];
            }
        }

        rzOld = rzNew;
    }

    // Step 4: Apply final displacements
    updateDrudePositions(state, particles, displacements);

    return true;
}

void DrudeTCGv2::computeFieldsAtParents(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& fields
) {
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};

    // Initialize
    for (auto& field : fields) {
        field[0] = field[1] = field[2] = 0.0;
    }

    // Compute field at each parent position from all fixed charges
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parent = state.atoms[particle.parentIndex];

        // Sum over all atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self and Drude
            if (j == particle.parentIndex || j == particle.drudeIndex) continue;

            const auto& atom = state.atoms[j];

            // Skip if no charge
            if (std::abs(atom.charge) < 1e-6) continue;

            // Skip other Drude particles
            bool isDrude = false;
            for (const auto& p : particles) {
                if (j == p.drudeIndex) {
                    isDrude = true;
                    break;
                }
            }
            if (isDrude) continue;

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

void DrudeTCGv2::applyMatrix(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& displacements,
    std::vector<Vec3>& result
) {
    // A*d = d - (q/k)*T*d
    // First term: d
    result = displacements;

    // Second term: (q/k)*T*d
    // T*d gives the induced fields due to the displaced Drudes
    std::vector<Vec3> inducedFields(particles.size());
    computeInducedFields(state, particles, inducedFields);

    // Subtract (q/k)*induced_fields
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        double factor = particle.charge / particle.kSpring;
        result[i][0] -= factor * inducedFields[i][0];
        result[i][1] -= factor * inducedFields[i][1];
        result[i][2] -= factor * inducedFields[i][2];
    }
}

void DrudeTCGv2::computeInducedFields(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& fields
) {
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};

    // Initialize
    for (auto& field : fields) {
        field[0] = field[1] = field[2] = 0.0;
    }

    // Compute field at each Drude position from other Drude charges
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle_i = particles[i];
        const auto& drude_i = state.atoms[particle_i.drudeIndex];

        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;

            const auto& particle_j = particles[j];
            const auto& drude_j = state.atoms[particle_j.drudeIndex];

            // Distance with PBC
            double dx = drude_i.x - drude_j.x;
            double dy = drude_i.y - drude_j.y;
            double dz = drude_i.z - drude_j.z;

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

void DrudeTCGv2::updateDrudePositions(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<Vec3>& displacements
) {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parent = state.atoms[particle.parentIndex];
        auto& drude = state.atoms[particle.drudeIndex];

        drude.x = parent.x + displacements[i][0];
        drude.y = parent.y + displacements[i][1];
        drude.z = parent.z + displacements[i][2];
    }
}

void DrudeTCGv2::buildPreconditioner(
    const std::vector<DrudeParticle>& particles
) {
    // Simple diagonal preconditioner: M_ii = 1 / (1 + alpha_i)
    // where alpha_i is a rough estimate of diagonal dominance
    for (size_t i = 0; i < particles.size(); ++i) {
        // The diagonal element of A is 1
        // The off-diagonal contributions scale with q/k
        // For stability, we use a conservative estimate
        precond_[i] = 1.0;  // Can be tuned based on system
    }
}

void DrudeTCGv2::applyPreconditioner(
    const std::vector<Vec3>& residuals,
    std::vector<Vec3>& z
) {
    for (size_t i = 0; i < residuals.size(); ++i) {
        z[i][0] = precond_[i] * residuals[i][0];
        z[i][1] = precond_[i] * residuals[i][1];
        z[i][2] = precond_[i] * residuals[i][2];
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
