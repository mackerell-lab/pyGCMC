#include "DrudeExperimentalCore.hpp"
#include "DrudeSCFOM.hpp"
#include "DrudeSCFOpenMMExact.hpp"
#include "DrudeNBTholeBuilder.hpp"
#include <algorithm>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

DrudeExperimentalCore::DrudeExperimentalCore() 
    : m_includeCoulomb(false) {
    m_scf = std::make_unique<DrudeSCFOM>();
    
    // Default to S1_POINT_CHARGE algorithm - the standard CHARMM/OpenMM model
    m_scf->setAlgorithm(DrudeAlgorithm::S1_POINT_CHARGE);
    
    // Set parameters optimized for standard S1 point charge model
    m_params.tolerance = 1e-5;
    m_params.maxIterations = 300;  // More iterations for better convergence
    m_params.dampingFactor = 0.9;  // Higher value = less damping, better convergence
    m_params.enableHardWall = false;  // Production default
    m_params.maxDrudeDistance = 0.02; // 0.2 Å
    m_params.excludePartnerParentInExternalField = false;  // Standard setting
}

double DrudeExperimentalCore::calculateEnergy(model::MCState& state) {
    if (m_particles.empty()) {
        return 0.0;
    }

    // Route to appropriate optimizer based on requireExactMatch flag
    bool converged = false;
    if (m_params.requireExactMatch) {
        // Use OpenMM-exact optimizer with dS/dr term
        DrudeSCFOpenMMExact exact;
        DrudeSCFOpenMMExactParams exactParams;
        exactParams.minimizationErrorTolerance = m_params.tolerance;
        exactParams.maxDrudeDistance = m_params.maxDrudeDistance;
        exactParams.maxIterations = m_params.maxIterations;
        exactParams.logLevel = m_params.logLevel;
        converged = exact.optimize(state, m_particles, m_screenedPairs, exactParams);
    } else {
        // Use fast SCF optimizer
        converged = m_scf->optimize(state, m_particles, m_screenedPairs, m_params);
    }
    
    if (!converged) {
        if (m_params.requireConvergence) {
            throw std::runtime_error("Experimental Drude SCF did not converge in calculateEnergy");
        }
        if (m_params.logLevel > 0) {
            std::cerr << "Warning: Experimental Drude SCF did not converge in calculateEnergy\n";
        }
    }

    // Return polarization (spring) energy only to avoid double counting
    double energy = 0.0;
    
    for (const auto& p : m_particles) {
        if (p.drudeIndex < 0 || p.drudeIndex >= state.activeAtomCount ||
            p.parentIndex < 0 || p.parentIndex >= state.activeAtomCount) {
            continue;
        }
        
        const auto& d = state.atoms[p.drudeIndex];
        const auto& o = state.atoms[p.parentIndex];
        
        double dx = d.x - o.x;
        double dy = d.y - o.y;
        double dz = d.z - o.z;
        
        // Apply PBC
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        DrudeSCFOM::applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        energy += 0.5 * p.kSpring * r2;
        
        // Add anisotropic terms if present
        if (p.aniso12 != 0 || p.aniso34 != 0) {
            // Simple anisotropy for now (can be expanded)
            // This would need proper orientation vectors in full implementation
        }
    }
    
    // Optionally include Coulomb energy for testing (normally should be false)
    if (m_includeCoulomb) {
        // This would add Coulomb energy calculation
        // Not implemented here as we want to avoid double counting
    }
    
    return energy;
}

void DrudeExperimentalCore::calculateForces(model::MCState& state, std::vector<Vec3>& forces) {
    if (forces.size() != static_cast<size_t>(state.activeAtomCount)) {
        forces.assign(state.activeAtomCount, {0.0, 0.0, 0.0});
    }
    
    if (m_particles.empty()) {
        return;
    }
    
    // Ensure positions are at SCF minimum
    bool converged = m_scf->optimize(state, m_particles, m_screenedPairs, m_params);
    if (!converged) {
        if (m_params.requireConvergence) {
            throw std::runtime_error("Experimental Drude SCF did not converge in calculateForces");
        }
        if (m_params.logLevel > 0) {
            std::cerr << "Warning: Experimental Drude SCF did not converge in calculateForces\n";
        }
    }
    
    // Spring forces only; Coulomb handled by main nonbonded module
    for (const auto& p : m_particles) {
        if (p.drudeIndex < 0 || p.drudeIndex >= state.activeAtomCount ||
            p.parentIndex < 0 || p.parentIndex >= state.activeAtomCount) {
            continue;
        }
        
        const auto& d = state.atoms[p.drudeIndex];
        const auto& o = state.atoms[p.parentIndex];
        
        double dx = d.x - o.x;
        double dy = d.y - o.y;
        double dz = d.z - o.z;
        
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        DrudeSCFOM::applyPBC(dx, dy, dz, box);
        
        Vec3 f = {-p.kSpring * dx, -p.kSpring * dy, -p.kSpring * dz};
        
        // Apply equal and opposite forces
        forces[p.drudeIndex][0] += f[0];
        forces[p.drudeIndex][1] += f[1];
        forces[p.drudeIndex][2] += f[2];
        
        forces[p.parentIndex][0] -= f[0];
        forces[p.parentIndex][1] -= f[1];
        forces[p.parentIndex][2] -= f[2];
    }
}

int DrudeExperimentalCore::addParticle(const DrudeParticle& particle) {
    auto p = particle;
    p.computeSpringConstants();
    m_particles.push_back(p);
    return static_cast<int>(m_particles.size()) - 1;
}

void DrudeExperimentalCore::addScreenedPair(const ScreenedPair& pair) {
    if (pair.dipole1 < 0 || pair.dipole2 < 0 || 
        static_cast<size_t>(pair.dipole1) >= m_particles.size() || 
        static_cast<size_t>(pair.dipole2) >= m_particles.size()) {
        throw std::runtime_error("Invalid screened pair indices");
    }
    m_screenedPairs.push_back(pair);
}

void DrudeExperimentalCore::autoScreenPairs(double thole, double cutoff_nm) {
    buildNBTholePairs(m_particles, thole, cutoff_nm, m_screenedPairs);
}

void DrudeExperimentalCore::setAlgorithm(platform::cpu::DrudeAlgorithm /*algorithm*/) {
    // Fixed to OpenMM-style SCF for experimental version
}

void DrudeExperimentalCore::setDrudeAlgorithm(int algo) {
    // Set the experimental algorithm
    // 0 = S1_POINT_CHARGE (CHARMM standard)
    // 1 = S3S5_DIPOLE_TENSOR
    // 2 = S1_DIPOLE_FIELD (hybrid)
    // 3 = DIRECT_COULOMB (no screening)
    exp::DrudeAlgorithm expAlgo;
    switch (algo) {
        case 0:
            expAlgo = exp::DrudeAlgorithm::S1_POINT_CHARGE;
            break;
        case 1:
            // DEPRECATED: S3S5_DIPOLE_TENSOR - use standard S1 instead
            expAlgo = exp::DrudeAlgorithm::S1_POINT_CHARGE;
            break;
        case 2:
            // DEPRECATED: S1_DIPOLE_FIELD - use standard S1 instead
            expAlgo = exp::DrudeAlgorithm::S1_POINT_CHARGE;
            break;
        case 3:
            expAlgo = exp::DrudeAlgorithm::DIRECT_COULOMB;
            break;
        default:
            expAlgo = exp::DrudeAlgorithm::S1_POINT_CHARGE;
            break;
    }
    m_scf->setAlgorithm(expAlgo);
}

int DrudeExperimentalCore::getDrudeAlgorithm() const {
    exp::DrudeAlgorithm expAlgo = m_scf->getAlgorithm();
    switch (expAlgo) {
        case exp::DrudeAlgorithm::S1_POINT_CHARGE:
            return 0;
        /* DEPRECATED: S3S5 and S1_DIPOLE_FIELD no longer exist
        case exp::DrudeAlgorithm::S3S5_DIPOLE_TENSOR:
            return 1;
        case exp::DrudeAlgorithm::S1_DIPOLE_FIELD:
            return 2;
        */
        case exp::DrudeAlgorithm::DIRECT_COULOMB:
            return 3;
        default:
            return 0;
    }
}

void DrudeExperimentalCore::setParameters(const DrudeSCFParams& params) {
    m_params = params;
}

void DrudeExperimentalCore::clear() {
    m_particles.clear();
    m_screenedPairs.clear();
}

bool DrudeExperimentalCore::inSameMolecule(int atom1, int atom2, const model::MCState& state) const {
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
            return true;  // Both atoms in same residue
        }
    }
    
    return false;  // Atoms in different residues
}

double DrudeExperimentalCore::getSpectralRadius(const model::MCState& state) const {
    // Get spectral radius from the SCF optimizer
    if (m_scf) {
        return m_scf->estimateSpectralRadius(state, m_particles, m_screenedPairs);
    }
    return 0.0;
}

int DrudeExperimentalCore::getSCFIterationCount() const {
    return m_scf ? m_scf->getIterationCount() : 0;
}

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc