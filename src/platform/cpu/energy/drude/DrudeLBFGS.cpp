/**
 * @file DrudeLBFGS.cpp
 * @brief L-BFGS optimizer implementation for Drude oscillators
 */

#include "DrudeLBFGS.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

// StateVector implementation
double DrudeLBFGS::StateVector::dot(const StateVector& other) const {
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        result += x[i] * other.x[i];
    }
    return result;
}

void DrudeLBFGS::StateVector::axpy(double a, const StateVector& y) {
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += a * y.x[i];
    }
}

void DrudeLBFGS::StateVector::scale(double a) {
    for (auto& xi : x) {
        xi *= a;
    }
}

double DrudeLBFGS::StateVector::norm() const {
    return std::sqrt(dot(*this));
}

// Main optimization function
bool DrudeLBFGS::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    if (particles.empty()) return true;

    // Initialize state vector
    StateVector x(particles.size());
    packState(state, particles, x);

    // Working vectors
    StateVector gradient(particles.size());
    StateVector gradient_prev(particles.size());
    StateVector direction(particles.size());

    // L-BFGS history
    std::deque<HistoryEntry> history;

    // Initial evaluation
    double energy = evaluateEnergy(x, particles, screenedPairs, state);
    evaluateGradient(x, particles, screenedPairs, state, gradient);

    // Check initial convergence
    double grad_norm = gradient.norm();
    if (grad_norm < params.tolerance) {
        unpackState(x, particles, state);
        return true;
    }

    // L-BFGS iterations
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        // Store previous gradient
        gradient_prev = gradient;

        // Compute search direction
        double H0 = 1.0;  // Initial Hessian approximation
        if (!history.empty()) {
            // Estimate H0 from most recent update
            const auto& recent = history.back();
            H0 = recent.s.dot(recent.y) / recent.y.dot(recent.y);
        }

        computeSearchDirection(gradient, history, H0, direction);

        // Line search
        StateVector x_prev = x;
        double step_size = lineSearch(x, direction, energy, gradient,
                                     particles, screenedPairs, state, gradient);

        if (step_size < 1e-10) {
            // Line search failed - try steepest descent
            direction = gradient;
            direction.scale(-1.0);
            x = x_prev;
            step_size = lineSearch(x, direction, energy, gradient,
                                  particles, screenedPairs, state, gradient);

            if (step_size < 1e-10) {
                // Cannot make progress
                unpackState(x, particles, state);
                return false;
            }
        }

        // Update energy
        energy = evaluateEnergy(x, particles, screenedPairs, state);

        // Update history
        HistoryEntry entry;
        entry.s = x;
        entry.s.axpy(-1.0, x_prev);  // s = x - x_prev

        entry.y = gradient;
        entry.y.axpy(-1.0, gradient_prev);  // y = g - g_prev

        double ys = entry.y.dot(entry.s);
        if (ys > 1e-10) {  // Ensure positive definiteness
            entry.rho = 1.0 / ys;
            history.push_back(std::move(entry));

            // Limit history size
            if (history.size() > MEMORY_SIZE) {
                history.pop_front();
            }
        }

        // Check convergence
        grad_norm = gradient.norm();
        if (grad_norm < params.tolerance) {
            unpackState(x, particles, state);
            return true;
        }
    }

    // Did not converge
    unpackState(x, particles, state);
    return false;
}

void DrudeLBFGS::packState(const model::MCState& state,
                           const std::vector<DrudeParticle>& particles,
                           StateVector& x) const {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& drude = state.atoms[particles[i].drudeIndex];
        const auto& parent = state.atoms[particles[i].parentIndex];

        // Store relative positions
        x[3*i + 0] = drude.x - parent.x;
        x[3*i + 1] = drude.y - parent.y;
        x[3*i + 2] = drude.z - parent.z;

        // Apply PBC
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(x[3*i], x[3*i+1], x[3*i+2], box);
    }
}

void DrudeLBFGS::unpackState(const StateVector& x,
                             const std::vector<DrudeParticle>& particles,
                             model::MCState& state) const {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& parent = state.atoms[particles[i].parentIndex];
        auto& drude = state.atoms[particles[i].drudeIndex];

        // Update absolute positions
        drude.x = parent.x + x[3*i + 0];
        drude.y = parent.y + x[3*i + 1];
        drude.z = parent.z + x[3*i + 2];
    }
}

double DrudeLBFGS::evaluateEnergy(const StateVector& x,
                                  const std::vector<DrudeParticle>& particles,
                                  const std::vector<ScreenedPair>& screenedPairs,
                                  model::MCState& state) const {
    // Update positions
    unpackState(x, particles, state);

    double energy = 0.0;

    // 1. Spring energy
    for (size_t i = 0; i < particles.size(); ++i) {
        double dx = x[3*i + 0];
        double dy = x[3*i + 1];
        double dz = x[3*i + 2];
        double r2 = dx*dx + dy*dy + dz*dz;

        energy += 0.5 * particles[i].kSpring * r2;
    }

    // 2. Drude-external charge interactions
    std::vector<Vec3> electricField(particles.size(), {0.0, 0.0, 0.0});
    calculateExternalField(state, particles, electricField);

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        // Potential energy = -q * E * d
        energy -= particle.charge * (
            electricField[i][0] * x[3*i + 0] +
            electricField[i][1] * x[3*i + 1] +
            electricField[i][2] * x[3*i + 2]
        );
    }

    // 3. Drude-Drude interactions (with Thole screening)
    for (const auto& pair : screenedPairs) {
        const auto& p1 = particles[pair.dipole1];
        const auto& p2 = particles[pair.dipole2];

        // Get parent positions
        const auto& parent1 = state.atoms[p1.parentIndex];
        const auto& parent2 = state.atoms[p2.parentIndex];

        // Calculate distance between parents
        double dx = parent1.x - parent2.x;
        double dy = parent1.y - parent2.y;
        double dz = parent1.z - parent2.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);

        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-12) continue;

        // Calculate dipole moments
        Vec3 mu1 = {x[3*pair.dipole1 + 0] * p1.charge,
                    x[3*pair.dipole1 + 1] * p1.charge,
                    x[3*pair.dipole1 + 2] * p1.charge};
        Vec3 mu2 = {x[3*pair.dipole2 + 0] * p2.charge,
                    x[3*pair.dipole2 + 1] * p2.charge,
                    x[3*pair.dipole2 + 2] * p2.charge};

        // Thole screening
        double screening = computeTholeScreening(r, p1.polarizability,
                                               p2.polarizability, pair.thole);

        // Dipole-dipole interaction energy
        double r3 = r * r * r;
        double mu1_dot_mu2 = mu1[0]*mu2[0] + mu1[1]*mu2[1] + mu1[2]*mu2[2];
        double mu1_dot_r = (mu1[0]*dx + mu1[1]*dy + mu1[2]*dz) / r;
        double mu2_dot_r = (mu2[0]*dx + mu2[1]*dy + mu2[2]*dz) / r;

        energy += DrudeConstants::ONE_4PI_EPS0 * screening * (
            mu1_dot_mu2 / r3 - 3.0 * mu1_dot_r * mu2_dot_r / (r3 * r*r)
        );
    }

    return energy;
}

void DrudeLBFGS::evaluateGradient(const StateVector& x,
                                  const std::vector<DrudeParticle>& particles,
                                  const std::vector<ScreenedPair>& screenedPairs,
                                  model::MCState& state,
                                  StateVector& gradient) const {
    // Update positions
    unpackState(x, particles, state);

    // Calculate electric field
    std::vector<Vec3> electricField(particles.size(), {0.0, 0.0, 0.0});
    calculateElectricField(state, particles, screenedPairs, electricField);

    // Calculate gradient = -Force
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];

        // Spring force: F = -k * displacement
        double fx = -particle.kSpring * x[3*i + 0];
        double fy = -particle.kSpring * x[3*i + 1];
        double fz = -particle.kSpring * x[3*i + 2];

        // Electric force: F = q * E
        fx += particle.charge * electricField[i][0];
        fy += particle.charge * electricField[i][1];
        fz += particle.charge * electricField[i][2];

        // Gradient is negative force
        gradient[3*i + 0] = -fx;
        gradient[3*i + 1] = -fy;
        gradient[3*i + 2] = -fz;
    }
}

void DrudeLBFGS::computeSearchDirection(const StateVector& gradient,
                                       const std::deque<HistoryEntry>& history,
                                       double H0,
                                       StateVector& direction) const {
    // L-BFGS two-loop recursion
    direction = gradient;

    // First loop
    std::vector<double> alpha(history.size());

    for (int i = history.size() - 1; i >= 0; --i) {
        const auto& entry = history[i];
        alpha[i] = entry.rho * entry.s.dot(direction);
        direction.axpy(-alpha[i], entry.y);
    }

    // Scale by initial Hessian
    direction.scale(H0);

    // Second loop
    for (size_t i = 0; i < history.size(); ++i) {
        const auto& entry = history[i];
        double beta = entry.rho * entry.y.dot(direction);
        direction.axpy(alpha[i] - beta, entry.s);
    }

    // Negate for descent direction
    direction.scale(-1.0);
}

double DrudeLBFGS::lineSearch(StateVector& x,
                              const StateVector& direction,
                              double f0,
                              const StateVector& g0,
                              const std::vector<DrudeParticle>& particles,
                              const std::vector<ScreenedPair>& screenedPairs,
                              model::MCState& state,
                              StateVector& gradient) const {
    // Backtracking line search with Armijo condition
    double alpha = 1.0;
    double dg0 = direction.dot(g0);

    if (dg0 > 0) {
        // Not a descent direction
        return 0.0;
    }

    StateVector x_new = x;

    for (int iter = 0; iter < MAX_LINE_SEARCH; ++iter) {
        // Try step
        x_new = x;
        x_new.axpy(alpha, direction);

        // Evaluate function
        double f_new = evaluateEnergy(x_new, particles, screenedPairs, state);

        // Check Armijo condition
        if (f_new <= f0 + ARMIJO_C1 * alpha * dg0) {
            // Accept step
            x = x_new;
            evaluateGradient(x, particles, screenedPairs, state, gradient);
            return alpha;
        }

        // Reduce step size
        alpha *= 0.5;

        if (alpha < 1e-10) {
            break;
        }
    }

    return 0.0;  // Line search failed
}

// Electric field calculation methods (reuse from SCF implementation)
void DrudeLBFGS::calculateElectricField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    std::vector<Vec3>& electricField
) const {
    // Field from external charges
    calculateExternalField(state, particles, electricField);

    // Field from induced dipoles
    calculateInducedField(state, particles, screenedPairs, electricField);
}

// Copy implementation from DrudeSCF.cpp for these methods
void DrudeLBFGS::calculateExternalField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& electricField
) const {
    // Implementation identical to DrudeSCF::calculateExternalField
    // (Copy the full implementation here)
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& drudeAtom = state.atoms[particle.drudeIndex];

        Vec3 field = {0.0, 0.0, 0.0};

        for (int j = 0; j < state.activeAtomCount; ++j) {
            if (j == particle.drudeIndex) continue;
            if (j == particle.parentIndex) continue;
            if (inSameMolecule(particle.drudeIndex, j, state)) continue;

            const auto& atom = state.atoms[j];

            double dx = drudeAtom.x - atom.x;
            double dy = drudeAtom.y - atom.y;
            double dz = drudeAtom.z - atom.z;
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(dx, dy, dz, box);

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;

            double r = std::sqrt(r2);
            double r3 = r2 * r;

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

void DrudeLBFGS::calculateInducedField(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    std::vector<Vec3>& electricField
) const {
    // Field from Drude-Drude interactions with Thole screening
    for (const auto& pair : screenedPairs) {
        const auto& particle1 = particles[pair.dipole1];
        const auto& particle2 = particles[pair.dipole2];

        // Get Drude positions
        const auto& drude1 = state.atoms[particle1.drudeIndex];
        const auto& drude2 = state.atoms[particle2.drudeIndex];
        const auto& parent1 = state.atoms[particle1.parentIndex];
        const auto& parent2 = state.atoms[particle2.parentIndex];

        // Calculate dipole moments
        Vec3 mu1 = {
            (drude1.x - parent1.x) * particle1.charge,
            (drude1.y - parent1.y) * particle1.charge,
            (drude1.z - parent1.z) * particle1.charge
        };

        Vec3 mu2 = {
            (drude2.x - parent2.x) * particle2.charge,
            (drude2.y - parent2.y) * particle2.charge,
            (drude2.z - parent2.z) * particle2.charge
        };

        // Calculate distance between dipole centers (parent atoms)
        double dx = parent1.x - parent2.x;
        double dy = parent1.y - parent2.y;
        double dz = parent1.z - parent2.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);

        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-12) continue;

        double r = std::sqrt(r2);
        double r3 = r2 * r;
        double r5 = r3 * r2;

        // Thole screening
        double screening = computeTholeScreening(r, particle1.polarizability,
                                               particle2.polarizability, pair.thole);

        // Dipole field tensor T_ij
        double rhat_x = dx / r;
        double rhat_y = dy / r;
        double rhat_z = dz / r;

        // Field at dipole 1 due to dipole 2
        double mu_dot_r = mu2[0]*rhat_x + mu2[1]*rhat_y + mu2[2]*rhat_z;

        double factor1 = DrudeConstants::ONE_4PI_EPS0 * screening / r3;
        double factor2 = DrudeConstants::ONE_4PI_EPS0 * screening * 3.0 * mu_dot_r / r5;

        electricField[pair.dipole1][0] += factor1 * mu2[0] - factor2 * dx;
        electricField[pair.dipole1][1] += factor1 * mu2[1] - factor2 * dy;
        electricField[pair.dipole1][2] += factor1 * mu2[2] - factor2 * dz;

        // Field at dipole 2 due to dipole 1 (opposite direction)
        mu_dot_r = mu1[0]*(-rhat_x) + mu1[1]*(-rhat_y) + mu1[2]*(-rhat_z);
        factor2 = DrudeConstants::ONE_4PI_EPS0 * screening * 3.0 * mu_dot_r / r5;

        electricField[pair.dipole2][0] += factor1 * mu1[0] + factor2 * dx;
        electricField[pair.dipole2][1] += factor1 * mu1[1] + factor2 * dy;
        electricField[pair.dipole2][2] += factor1 * mu1[2] + factor2 * dz;
    }
}

void DrudeLBFGS::applyPBC(double& dx, double& dy, double& dz,
                          const std::array<double, 3>& box) const {
    dx -= box[0] * std::round(dx / box[0]);
    dy -= box[1] * std::round(dy / box[1]);
    dz -= box[2] * std::round(dz / box[2]);
}

bool DrudeLBFGS::inSameMolecule(int atom1, int atom2,
                                const model::MCState& state) const {
    // Simple check: atoms in same residue are in same molecule
    for (int i = 0; i < state.activeResidueCount; ++i) {
        const auto& res = state.residues[i];
        bool atom1InRes = (atom1 >= res.atomStart &&
                          atom1 < res.atomStart + res.atomCount);
        bool atom2InRes = (atom2 >= res.atomStart &&
                          atom2 < res.atomStart + res.atomCount);
        if (atom1InRes && atom2InRes) {
            return true;
        }
    }
    return false;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
