/**
 * @file DrudeSCFOpenMMExact.cpp
 * @brief OpenMM-exact Drude SCF implementation
 */

#include "DrudeSCFOpenMMExact.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <map>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

bool DrudeSCFOpenMMExact::optimize(model::montecarlo::MCState& state,
                                   const std::vector<DrudeParticle>& particles,
                                   const std::vector<ScreenedPair>& pairs,
                                   const DrudeSCFOpenMMExactParams& params) {
    
    if (particles.empty()) {
        m_lastIterationCount = 0;
        return true;
    }
    
    if (params.logLevel >= 1) {
        std::cout << "DrudeSCFOpenMMExact: Starting optimization with "
                  << particles.size() << " Drude particles" << std::endl;
    }
    
    // Initialize convergence tracking
    m_lastIterationCount = 0;
    double prevEnergy = calculateTotalEnergy(state, particles, pairs);
    m_lastEnergy = prevEnergy;
    
    // Force vector
    std::vector<Vec3> forces(particles.size());
    
    // Optimization loop (gradient descent with line search)
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        m_lastIterationCount = iter + 1;
        
        // Calculate forces (negative gradient)
        calculateForces(state, particles, pairs, forces);
        
        // Calculate force norm for convergence check
        double forceNorm = 0.0;
        for (const auto& f : forces) {
            forceNorm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
        }
        forceNorm = std::sqrt(forceNorm);
        m_lastForceNorm = forceNorm;
        
        if (params.logLevel >= 2) {
            std::cout << "  Iter " << iter << ": E=" << std::scientific 
                     << std::setprecision(8) << prevEnergy 
                     << ", |F|=" << forceNorm << std::endl;
        }
        
        // Check force convergence (OpenMM criterion)
        if (forceNorm < params.minimizationErrorTolerance) {
            if (params.logLevel >= 1) {
                std::cout << "Converged in " << iter + 1 
                         << " iterations (|F|=" << forceNorm << ")" << std::endl;
            }
            return true;
        }
        
        // Line search for optimal step
        double stepSize = lineSearch(state, particles, pairs, forces, params.stepSize);
        
        // Update positions
        if (!updatePositions(state, particles, forces, stepSize, params)) {
            if (params.logLevel >= 1) {
                std::cout << "Warning: Position update failed at iteration " 
                         << iter << std::endl;
            }
            return false;
        }
        
        // Calculate new energy
        double newEnergy = calculateTotalEnergy(state, particles, pairs);
        
        // Check energy convergence
        double energyChange = std::abs(newEnergy - prevEnergy);
        if (energyChange < params.energyTolerance && forceNorm < params.minimizationErrorTolerance * 10) {
            if (params.logLevel >= 1) {
                std::cout << "Converged by energy criterion in " << iter + 1 
                         << " iterations" << std::endl;
            }
            m_lastEnergy = newEnergy;
            return true;
        }
        
        prevEnergy = newEnergy;
        m_lastEnergy = newEnergy;
    }
    
    if (params.logLevel >= 1) {
        std::cout << "Warning: Did not converge in " << params.maxIterations 
                 << " iterations (|F|=" << m_lastForceNorm << ")" << std::endl;
    }
    
    return false;
}

double DrudeSCFOpenMMExact::calculateTotalEnergy(const model::montecarlo::MCState& state,
                                                 const std::vector<DrudeParticle>& particles,
                                                 const std::vector<ScreenedPair>& pairs) const {
    double energy = 0.0;
    
    // 1. Spring energy
    energy += calculateSpringEnergy(state, particles);
    
    // 2. Unscreened Coulomb energy (external field)
    energy += calculateUnscreenedEnergy(state, particles);
    
    // 3. Screened dipole-dipole interactions
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    
    // Build screened pairs map
    std::map<std::pair<int,int>, double> screenedPairsMap;
    for (const auto& pair : pairs) {
        screenedPairsMap[{pair.dipole1, pair.dipole2}] = pair.thole;
        screenedPairsMap[{pair.dipole2, pair.dipole1}] = pair.thole;
    }
    
    // Calculate all dipole-dipole interactions
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            const auto& pi = particles[i];
            const auto& pj = particles[j];
            
            // Get Thole parameter
            double thole = 0.0;
            auto it = screenedPairsMap.find({static_cast<int>(i), static_cast<int>(j)});
            if (it != screenedPairsMap.end()) {
                thole = it->second;
            }
            
            // Skip if no screening
            if (thole < 1e-10) continue;
            
            const auto& parent1 = state.atoms[pi.parentIndex];
            const auto& drude1 = state.atoms[pi.drudeIndex];
            const auto& parent2 = state.atoms[pj.parentIndex];
            const auto& drude2 = state.atoms[pj.drudeIndex];
            
            // P1-P2 (no screening in OpenMM)
            // Handled in unscreened energy
            
            // P1-D2 (screened)
            {
                Vec3 r1 = {parent1.x, parent1.y, parent1.z};
                Vec3 r2 = {drude2.x, drude2.y, drude2.z};
                double alpha_eff = std::cbrt(pj.polarizability);  // α^(1/3) for P-D
                energy += calculateScreenedEnergy(r1, r2, parent1.charge, pj.charge,
                                                 alpha_eff, thole, box);
            }
            
            // D1-P2 (screened)
            {
                Vec3 r1 = {drude1.x, drude1.y, drude1.z};
                Vec3 r2 = {parent2.x, parent2.y, parent2.z};
                double alpha_eff = std::cbrt(pi.polarizability);  // α^(1/3) for P-D
                energy += calculateScreenedEnergy(r1, r2, pi.charge, parent2.charge,
                                                 alpha_eff, thole, box);
            }
            
            // D1-D2 (screened)
            {
                Vec3 r1 = {drude1.x, drude1.y, drude1.z};
                Vec3 r2 = {drude2.x, drude2.y, drude2.z};
                double alpha_eff = std::pow(pi.polarizability * pj.polarizability, 1.0/6.0);
                energy += calculateScreenedEnergy(r1, r2, pi.charge, pj.charge,
                                                 alpha_eff, thole, box);
            }
        }
    }
    
    return energy;
}

void DrudeSCFOpenMMExact::calculateForces(const model::montecarlo::MCState& state,
                                          const std::vector<DrudeParticle>& particles,
                                          const std::vector<ScreenedPair>& pairs,
                                          std::vector<Vec3>& forces) const {
    // Clear forces
    for (auto& f : forces) {
        f[0] = f[1] = f[2] = 0.0;
    }
    
    // 1. Spring forces
    calculateSpringForces(state, particles, forces);
    
    // 2. Build necessary maps
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    std::vector<bool> isDrude(state.activeAtomCount, false);
    std::vector<int> atomToDipole(state.activeAtomCount, -1);
    
    for (size_t i = 0; i < particles.size(); ++i) {
        isDrude[particles[i].drudeIndex] = true;
        atomToDipole[particles[i].parentIndex] = i;
        atomToDipole[particles[i].drudeIndex] = i;
    }
    
    std::map<std::pair<int,int>, double> screenedPairsMap;
    for (const auto& pair : pairs) {
        screenedPairsMap[{pair.dipole1, pair.dipole2}] = pair.thole;
        screenedPairsMap[{pair.dipole2, pair.dipole1}] = pair.thole;
    }
    
    // 3. Calculate forces on each Drude particle
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& pi = particles[i];
        const auto& drude_i = state.atoms[pi.drudeIndex];
        Vec3 r_i = {drude_i.x, drude_i.y, drude_i.z};
        
        // External field forces (from non-dipole charges)
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self, parent, and other Drude particles
            if (j == pi.drudeIndex || j == pi.parentIndex || isDrude[j]) continue;
            
            // Skip if part of another dipole (parent)
            if (atomToDipole[j] >= 0) continue;
            
            // Skip same molecule
            if (inSameMolecule(pi.parentIndex, j, state)) continue;
            
            const auto& atom_j = state.atoms[j];
            if (std::abs(atom_j.charge) < 1e-14) continue;
            
            Vec3 r_j = {atom_j.x, atom_j.y, atom_j.z};
            
            // Unscreened Coulomb force
            double dx = r_i[0] - r_j[0];
            double dy = r_i[1] - r_j[1];
            double dz = r_i[2] - r_j[2];
            applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > 1e-12) {
                double r = std::sqrt(r2);
                double r3 = r2 * r;
                double factor = DrudeSCFOpenMMExactParams::ONE_4PI_EPS0 * 
                               pi.charge * atom_j.charge / r3;
                
                forces[i][0] -= factor * dx;  // Force is negative gradient
                forces[i][1] -= factor * dy;
                forces[i][2] -= factor * dz;
            }
        }
        
        // Induced field forces (from other dipoles)
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;
            
            const auto& pj = particles[j];
            
            // Get Thole parameter
            double thole = 0.0;
            auto it = screenedPairsMap.find({static_cast<int>(i), static_cast<int>(j)});
            if (it != screenedPairsMap.end()) {
                thole = it->second;
            }
            
            const auto& parent_j = state.atoms[pj.parentIndex];
            const auto& drude_j = state.atoms[pj.drudeIndex];
            
            // P-D force (parent_j on drude_i) with dS/dr
            {
                Vec3 r_j = {parent_j.x, parent_j.y, parent_j.z};
                double alpha_eff = std::cbrt(pi.polarizability);
                
                Vec3 force_pd = calculateScreenedForce(r_i, r_j, pi.charge, parent_j.charge,
                                                       alpha_eff, thole, box);
                forces[i][0] -= force_pd[0];
                forces[i][1] -= force_pd[1];
                forces[i][2] -= force_pd[2];
            }
            
            // D-D force (drude_j on drude_i) with dS/dr
            {
                Vec3 r_j = {drude_j.x, drude_j.y, drude_j.z};
                double alpha_eff = std::pow(pi.polarizability * pj.polarizability, 1.0/6.0);
                
                Vec3 force_dd = calculateScreenedForce(r_i, r_j, pi.charge, pj.charge,
                                                       alpha_eff, thole, box);
                forces[i][0] -= force_dd[0];
                forces[i][1] -= force_dd[1];
                forces[i][2] -= force_dd[2];
            }
        }
    }
}

double DrudeSCFOpenMMExact::calculateScreenedEnergy(const Vec3& r_i, const Vec3& r_j,
                                                    double q_i, double q_j,
                                                    double alpha_eff, double thole,
                                                    const std::array<double, 3>& box) const {
    double dx = r_i[0] - r_j[0];
    double dy = r_i[1] - r_j[1];
    double dz = r_i[2] - r_j[2];
    
    applyPBC(dx, dy, dz, box);
    
    double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 < 1e-12) return 0.0;
    
    double r = std::sqrt(r2);
    
    // Calculate S1
    double u = thole * r / alpha_eff;
    double S1 = (u > 50.0) ? 1.0 : 1.0 - (1.0 + 0.5*u) * std::exp(-u);
    
    // Screened Coulomb energy
    return DrudeSCFOpenMMExactParams::ONE_4PI_EPS0 * q_i * q_j * S1 / r;
}

DrudeSCFOpenMMExact::Vec3 DrudeSCFOpenMMExact::calculateScreenedForce(const Vec3& r_i, const Vec3& r_j,
                                                 double q_i, double q_j,
                                                 double alpha_eff, double thole,
                                                 const std::array<double, 3>& box) const {
    Vec3 force = {0.0, 0.0, 0.0};
    
    double dx = r_i[0] - r_j[0];
    double dy = r_i[1] - r_j[1];
    double dz = r_i[2] - r_j[2];
    
    applyPBC(dx, dy, dz, box);
    
    double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 < 1e-12) return force;
    
    double r = std::sqrt(r2);
    double r3 = r2 * r;
    
    // Calculate u and S1 with derivative
    double u = thole * r / alpha_eff;
    S1Derivatives s1d = S1Derivatives::calculate(u);
    
    // du/dr
    double du_dr = thole / alpha_eff;
    
    // Complete force expression with dS/dr term
    // F = q_i * q_j * [S1/r³ - (dS1/du)(du/dr)/r²] * r_vec (OpenMM convention)
    double factor = DrudeSCFOpenMMExactParams::ONE_4PI_EPS0 * q_i * q_j;
    double term1 = s1d.S1 / r3;                    // S1/r³
    double term2 = -s1d.dS1_du * du_dr / r2;      // -(dS1/du)(du/dr)/r² [negative sign for OpenMM match]
    
    double total_factor = factor * (term1 + term2);
    
    force[0] = total_factor * dx;
    force[1] = total_factor * dy;
    force[2] = total_factor * dz;
    
    return force;
}

double DrudeSCFOpenMMExact::calculateSpringEnergy(const model::montecarlo::MCState& state,
                                                  const std::vector<DrudeParticle>& particles) const {
    double energy = 0.0;
    
    for (const auto& p : particles) {
        const auto& parent = state.atoms[p.parentIndex];
        const auto& drude = state.atoms[p.drudeIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        
        // Apply PBC for spring (same as OpenMM)
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        
        // Spring energy: U = 0.5 * k * r²
        energy += 0.5 * p.kSpring * r2;
    }
    
    return energy;
}

void DrudeSCFOpenMMExact::calculateSpringForces(const model::montecarlo::MCState& state,
                                                const std::vector<DrudeParticle>& particles,
                                                std::vector<Vec3>& forces) const {
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        const auto& parent = state.atoms[p.parentIndex];
        const auto& drude = state.atoms[p.drudeIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        
        applyPBC(dx, dy, dz, box);
        
        // Spring force: F = -k * r
        forces[i][0] -= p.kSpring * dx;
        forces[i][1] -= p.kSpring * dy;
        forces[i][2] -= p.kSpring * dz;
    }
}

double DrudeSCFOpenMMExact::calculateUnscreenedEnergy(const model::montecarlo::MCState& state,
                                                      const std::vector<DrudeParticle>& particles) const {
    double energy = 0.0;
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    
    // Mark Drude atoms
    std::vector<bool> isDrude(state.activeAtomCount, false);
    std::vector<int> atomToDipole(state.activeAtomCount, -1);
    
    for (size_t i = 0; i < particles.size(); ++i) {
        isDrude[particles[i].drudeIndex] = true;
        atomToDipole[particles[i].parentIndex] = i;
        atomToDipole[particles[i].drudeIndex] = i;
    }
    
    // Calculate unscreened Coulomb interactions
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& pi = particles[i];
        const auto& drude_i = state.atoms[pi.drudeIndex];
        
        // Interactions with external charges
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self, parent, and other Drude/parent atoms
            if (j == pi.drudeIndex || j == pi.parentIndex) continue;
            if (isDrude[j] || atomToDipole[j] >= 0) continue;
            
            // Skip same molecule
            if (inSameMolecule(pi.parentIndex, j, state)) continue;
            
            const auto& atom_j = state.atoms[j];
            if (std::abs(atom_j.charge) < 1e-14) continue;
            
            double dx = drude_i.x - atom_j.x;
            double dy = drude_i.y - atom_j.y;
            double dz = drude_i.z - atom_j.z;
            
            applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > 1e-12) {
                double r = std::sqrt(r2);
                energy += DrudeSCFOpenMMExactParams::ONE_4PI_EPS0 * 
                         pi.charge * atom_j.charge / r;
            }
        }
    }
    
    return energy;
}

bool DrudeSCFOpenMMExact::updatePositions(model::montecarlo::MCState& state,
                                          const std::vector<DrudeParticle>& particles,
                                          const std::vector<Vec3>& forces,
                                          double stepSize,
                                          const DrudeSCFOpenMMExactParams& params) const {
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        auto& drude = state.atoms[p.drudeIndex];
        const auto& parent = state.atoms[p.parentIndex];
        
        // Gradient descent step (force points in direction of decreasing energy)
        double new_x = drude.x + stepSize * forces[i][0];
        double new_y = drude.y + stepSize * forces[i][1];
        double new_z = drude.z + stepSize * forces[i][2];
        
        // Check max distance constraint
        double dx = new_x - parent.x;
        double dy = new_y - parent.y;
        double dz = new_z - parent.z;
        
        applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > params.maxDrudeDistance * params.maxDrudeDistance) {
            // Apply hard wall constraint
            double r = std::sqrt(r2);
            double scale = params.maxDrudeDistance / r;
            dx *= scale;
            dy *= scale;
            dz *= scale;
            
            new_x = parent.x + dx;
            new_y = parent.y + dy;
            new_z = parent.z + dz;
        }
        
        // Update position
        drude.x = new_x;
        drude.y = new_y;
        drude.z = new_z;
    }
    
    return true;
}

double DrudeSCFOpenMMExact::lineSearch(const model::montecarlo::MCState& state,
                                       const std::vector<DrudeParticle>& particles,
                                       const std::vector<ScreenedPair>& pairs,
                                       const std::vector<Vec3>& searchDirection,
                                       double initialStep) const {
    // Simple backtracking line search
    double alpha = initialStep;
    const double c = 0.5;  // Armijo constant
    const double rho = 0.5; // Backtracking factor
    
    // Calculate initial energy
    double E0 = calculateTotalEnergy(state, particles, pairs);
    
    // Calculate directional derivative (gradient · search direction)
    std::vector<Vec3> gradient(particles.size());
    calculateForces(state, particles, pairs, gradient);
    
    double dirDeriv = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        dirDeriv += gradient[i][0] * searchDirection[i][0] +
                   gradient[i][1] * searchDirection[i][1] +
                   gradient[i][2] * searchDirection[i][2];
    }
    
    // Save original positions
    std::vector<Vec3> originalPos;
    for (const auto& p : particles) {
        const auto& drude = state.atoms[p.drudeIndex];
        originalPos.push_back({drude.x, drude.y, drude.z});
    }
    
    // Backtracking loop
    const int maxBacktrack = 10;
    for (int k = 0; k < maxBacktrack; ++k) {
        // Try step with current alpha
        model::montecarlo::MCState testState = state;  // Copy state
        
        // Update positions with current step size
        for (size_t i = 0; i < particles.size(); ++i) {
            auto& drude = testState.atoms[particles[i].drudeIndex];
            drude.x = originalPos[i][0] + alpha * searchDirection[i][0];
            drude.y = originalPos[i][1] + alpha * searchDirection[i][1];
            drude.z = originalPos[i][2] + alpha * searchDirection[i][2];
        }
        
        // Calculate new energy
        double E_new = calculateTotalEnergy(testState, particles, pairs);
        
        // Check Armijo condition
        if (E_new <= E0 + c * alpha * dirDeriv) {
            return alpha;  // Accept this step size
        }
        
        // Reduce step size
        alpha *= rho;
    }
    
    // Return minimum step size if backtracking failed
    return alpha;
}

void DrudeSCFOpenMMExact::applyPBC(double& dx, double& dy, double& dz,
                                   const std::array<double, 3>& box) const {
    // Minimum image convention
    if (box[0] > 0) {
        dx -= box[0] * std::round(dx / box[0]);
    }
    if (box[1] > 0) {
        dy -= box[1] * std::round(dy / box[1]);
    }
    if (box[2] > 0) {
        dz -= box[2] * std::round(dz / box[2]);
    }
}

bool DrudeSCFOpenMMExact::inSameMolecule(int atom1, int atom2, 
                                         const model::montecarlo::MCState& state) const {
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
    
    return res1 == res2 && res1 >= 0;
}

double DrudeSCFOpenMMExact::calculateEnergy(const model::montecarlo::MCState& state,
                                           const std::vector<DrudeParticle>& particles,
                                           const std::vector<ScreenedPair>& pairs) const {
    return calculateTotalEnergy(state, particles, pairs);
}

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc