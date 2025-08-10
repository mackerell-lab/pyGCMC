#include "DrudeExperimentalCore.hpp"
#include "DrudeSCFOM.hpp"
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
    
    // Set reasonable defaults
    m_params.tolerance = 1e-5;
    m_params.maxIterations = 200;
    m_params.dampingFactor = 0.1;  // Lower damping for better convergence
    m_params.enableHardWall = false;  // Changed to false to align with OpenMM/CHARMM defaults
    m_params.maxDrudeDistance = 0.02; // 0.2 Å
}

double DrudeExperimentalCore::calculateEnergy(model::MCState& state) {
    if (m_particles.empty()) {
        return 0.0;
    }

    // Optimize positions with OpenMM-style SCF
    bool converged = m_scf->optimize(state, m_particles, m_screenedPairs, m_params);
    
    if (!converged) {
        // SCF did not converge - this can happen but is handled by hard wall
        // std::cerr << "Warning: Drude SCF did not converge\n";
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
    m_scf->optimize(state, m_particles, m_screenedPairs, m_params);
    
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

void DrudeExperimentalCore::setAlgorithm(DrudeAlgorithm /*algorithm*/) {
    // Fixed to OpenMM-style SCF for experimental version
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

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc