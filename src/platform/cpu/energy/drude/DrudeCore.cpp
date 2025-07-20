/**
 * @file DrudeCore.cpp
 * @brief Core implementation of Drude oscillator calculations
 */

#include "DrudeCore.hpp"
#include "DrudeStructures.hpp"
#include <cmath>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

// Global instance for static interface
static DrudeCore g_drudeCore;

DrudeCore::DrudeCore() 
    : m_algorithm(DrudeAlgorithm::SCF),
      m_scfOptimizer(std::make_unique<DrudeSCF>()),
      m_currentOptimizer(m_scfOptimizer.get()) {
}

double DrudeCore::calculateEnergy(model::MCState& state) {
    if (!m_currentOptimizer) {
        throw std::runtime_error("No Drude optimizer available");
    }
    
    // First optimize Drude positions
    bool converged = m_currentOptimizer->optimize(state, m_particles, m_screenedPairs, m_params);
    
    if (!converged) {
        std::cerr << "Warning: Drude SCF did not converge\n";
    }
    
    // Calculate total energy
    double energy = 0.0;
    
    // 1. Harmonic spring energy
    energy += calculateHarmonicEnergy(state);
    
    // 2. Thole-screened Coulomb energy
    energy += calculateScreenedCoulombEnergy(state);
    
    return energy;
}

void DrudeCore::calculateForces(model::MCState& state, std::vector<Vec3>& forces) {
    // Ensure forces vector has correct size
    if (forces.size() != static_cast<size_t>(state.activeAtomCount)) {
        forces.resize(static_cast<size_t>(state.activeAtomCount), {0.0, 0.0, 0.0});
    }
    
    // First optimize Drude positions
    m_currentOptimizer->optimize(state, m_particles, m_screenedPairs, m_params);
    
    // Calculate forces from harmonic springs
    calculateHarmonicForces(state, forces);
    
    // Calculate forces from Thole-screened interactions
    calculateScreenedCoulombForces(state, forces);
}

int DrudeCore::addParticle(const DrudeParticle& particle) {
    // Make a copy and compute derived quantities
    DrudeParticle p = particle;
    p.computeSpringConstants();
    
    m_particles.push_back(p);
    return m_particles.size() - 1;
}

void DrudeCore::addScreenedPair(const ScreenedPair& pair) {
    // Validate indices
    if (pair.dipole1 < 0 || static_cast<size_t>(pair.dipole1) >= m_particles.size() || 
        pair.dipole2 < 0 || static_cast<size_t>(pair.dipole2) >= m_particles.size()) {
        throw std::runtime_error("Invalid dipole indices in screened pair");
    }
    
    m_screenedPairs.push_back(pair);
}

void DrudeCore::setAlgorithm(DrudeAlgorithm algorithm) {
    if (algorithm == m_algorithm) return;
    
    m_algorithm = algorithm;
    
    switch (algorithm) {
        case DrudeAlgorithm::SCF:
            if (!m_scfOptimizer) {
                m_scfOptimizer = std::make_unique<DrudeSCF>();
            }
            m_currentOptimizer = m_scfOptimizer.get();
            break;
            
        case DrudeAlgorithm::OPT3:
            if (!m_opt3Optimizer) {
                m_opt3Optimizer = std::make_unique<DrudeOPT3>();
            }
            m_currentOptimizer = m_opt3Optimizer.get();
            break;
            
        case DrudeAlgorithm::FBP:
            if (!m_fbpOptimizer) {
                m_fbpOptimizer = std::make_unique<DrudeFBP>();
            }
            m_currentOptimizer = m_fbpOptimizer.get();
            break;
    }
}

void DrudeCore::setParameters(const DrudeSCFParams& params) {
    m_params = params;
}

void DrudeCore::clear() {
    m_particles.clear();
    m_screenedPairs.clear();
}

size_t DrudeCore::getNumParticles() const {
    return m_particles.size();
}

double DrudeCore::calculateHarmonicEnergy(const model::MCState& state) const {
    double energy = 0.0;
    
    for (const auto& particle : m_particles) {
        // Get positions
        const auto& drudeAtom = state.atoms[particle.drudeIndex];
        const auto& parentAtom = state.atoms[particle.parentIndex];
        
        // Displacement vector
        double dx = drudeAtom.x - parentAtom.x;
        double dy = drudeAtom.y - parentAtom.y;
        double dz = drudeAtom.z - parentAtom.z;
        
        // Apply periodic boundary conditions
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        // Isotropic harmonic energy: 0.5 * k * r^2
        double r2 = dx*dx + dy*dy + dz*dz;
        energy += 0.5 * particle.kSpring * r2;
        
        // Anisotropic contributions (if present)
        if (particle.aniso1Index >= 0 && particle.aniso2Index >= 0) {
            // Calculate projection along anisotropy axis
            const auto& aniso1 = state.atoms[particle.aniso1Index];
            const auto& aniso2 = state.atoms[particle.aniso2Index];
            
            double ax = aniso2.x - aniso1.x;
            double ay = aniso2.y - aniso1.y;
            double az = aniso2.z - aniso1.z;
            std::array<double, 3> box2 = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(ax, ay, az, box2);
            
            double anorm = std::sqrt(ax*ax + ay*ay + az*az);
            if (anorm > 1e-6) {
                ax /= anorm;
                ay /= anorm;
                az /= anorm;
                
                double proj = dx*ax + dy*ay + dz*az;
                energy += 0.5 * particle.kAniso1 * proj * proj;
            }
        }
        
        if (particle.aniso3Index >= 0 && particle.aniso4Index >= 0) {
            // Second anisotropy axis
            const auto& aniso3 = state.atoms[particle.aniso3Index];
            const auto& aniso4 = state.atoms[particle.aniso4Index];
            
            double ax = aniso4.x - aniso3.x;
            double ay = aniso4.y - aniso3.y;
            double az = aniso4.z - aniso3.z;
            std::array<double, 3> box2 = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(ax, ay, az, box2);
            
            double anorm = std::sqrt(ax*ax + ay*ay + az*az);
            if (anorm > 1e-6) {
                ax /= anorm;
                ay /= anorm;
                az /= anorm;
                
                double proj = dx*ax + dy*ay + dz*az;
                energy += 0.5 * particle.kAniso2 * proj * proj;
            }
        }
    }
    
    return energy;
}

double DrudeCore::calculateScreenedCoulombEnergy(const model::MCState& state) const {
    double energy = 0.0;
    
    for (const auto& pair : m_screenedPairs) {
        const auto& particle1 = m_particles[pair.dipole1];
        const auto& particle2 = m_particles[pair.dipole2];
        
        // Get Drude positions
        const auto& drude1 = state.atoms[particle1.drudeIndex];
        const auto& drude2 = state.atoms[particle2.drudeIndex];
        
        // Calculate Drude-Drude distance
        double dx = drude2.x - drude1.x;
        double dy = drude2.y - drude1.y;
        double dz = drude2.z - drude1.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-6) continue;  // Skip if too close
        
        // Calculate Thole screening
        double screening = computeTholeScreening(r, particle1.polarizability, 
                                                 particle2.polarizability, pair.thole);
        
        // Screened Coulomb energy
        double q1 = particle1.charge;
        double q2 = particle2.charge;
        energy += DrudeConstants::ONE_4PI_EPS0 * q1 * q2 * screening / r;
    }
    
    return energy;
}

void DrudeCore::calculateHarmonicForces(const model::MCState& state, 
                                        std::vector<Vec3>& forces) const {
    for (const auto& particle : m_particles) {
        // Get positions
        const auto& drudeAtom = state.atoms[particle.drudeIndex];
        const auto& parentAtom = state.atoms[particle.parentIndex];
        
        // Displacement vector
        double dx = drudeAtom.x - parentAtom.x;
        double dy = drudeAtom.y - parentAtom.y;
        double dz = drudeAtom.z - parentAtom.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        // Isotropic harmonic force: F = -k * r
        Vec3 force = {-particle.kSpring * dx,
                      -particle.kSpring * dy,
                      -particle.kSpring * dz};
        
        // Apply equal and opposite forces
        forces[particle.drudeIndex][0] += force[0];
        forces[particle.drudeIndex][1] += force[1];
        forces[particle.drudeIndex][2] += force[2];
        forces[particle.parentIndex][0] -= force[0];
        forces[particle.parentIndex][1] -= force[1];
        forces[particle.parentIndex][2] -= force[2];
        
        // Anisotropic forces (if present)
        // ... (similar to energy calculation but with force derivatives)
    }
}

void DrudeCore::calculateScreenedCoulombForces(const model::MCState& state,
                                               std::vector<Vec3>& forces) const {
    for (const auto& pair : m_screenedPairs) {
        const auto& particle1 = m_particles[pair.dipole1];
        const auto& particle2 = m_particles[pair.dipole2];
        
        // Get Drude positions
        const auto& drude1 = state.atoms[particle1.drudeIndex];
        const auto& drude2 = state.atoms[particle2.drudeIndex];
        
        // Calculate displacement
        double dx = drude2.x - drude1.x;
        double dy = drude2.y - drude1.y;
        double dz = drude2.z - drude1.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        double r = std::sqrt(r2);
        if (r < 1e-6) continue;
        
        // Calculate screening and its derivative
        double alpha_ij = std::pow(particle1.polarizability * particle2.polarizability, 1.0/6.0);
        double u = r / alpha_ij;
        double thu = pair.thole * u;
        double exp_thu = std::exp(-thu);
        
        double screening = 1.0 - (1.0 + thu/2.0) * exp_thu;
        double dscreening_dr = pair.thole * thu * exp_thu / (2.0 * alpha_ij);
        
        // Calculate force
        double q1q2 = particle1.charge * particle2.charge;
        double prefactor = DrudeConstants::ONE_4PI_EPS0 * q1q2;
        
        // F = -dE/dr * r_hat
        double f_mag = prefactor * (screening/r2 + dscreening_dr/r) / r;
        
        Vec3 force = {f_mag * dx, f_mag * dy, f_mag * dz};
        
        // Apply forces
        forces[particle1.drudeIndex][0] -= force[0];
        forces[particle1.drudeIndex][1] -= force[1];
        forces[particle1.drudeIndex][2] -= force[2];
        forces[particle2.drudeIndex][0] += force[0];
        forces[particle2.drudeIndex][1] += force[1];
        forces[particle2.drudeIndex][2] += force[2];
    }
}

void DrudeCore::applyPBC(double& dx, double& dy, double& dz, const std::array<double, 3>& box) const {
    dx -= box[0] * std::round(dx / box[0]);
    dy -= box[1] * std::round(dy / box[1]);
    dz -= box[2] * std::round(dz / box[2]);
}

// Static interface implementation
DrudeCore& DrudeCore::getInstance() {
    return g_drudeCore;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc