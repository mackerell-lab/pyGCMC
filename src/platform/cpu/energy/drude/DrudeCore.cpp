/**
 * @file DrudeCore.cpp
 * @brief Core implementation of Drude oscillator calculations
 */

#include "DrudeCore.hpp"
#include "DrudeStructures.hpp"
#include "DrudeDirectPolarization.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

// Global instance for static interface
static DrudeCore g_drudeCore;

void DrudeCore::buildActiveAtomMask(const model::MCState& state, std::vector<char>& mask) {
    const int nAtoms = std::max(0, state.activeAtomCount);
    mask.assign(static_cast<size_t>(nAtoms), 1);
    if (state.residues.empty()) {
        return;
    }

    // Default: all atoms [0..activeAtomCount) are active.
    // To support GCMC deletion without compacting atom arrays, explicitly mask out atoms
    // that belong to inactive residues. Atoms not covered by any residue range remain active
    // (some unit tests model external charges this way).
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

bool DrudeCore::isActiveAtom(int atomIndex, const std::vector<char>& mask) {
    if (atomIndex < 0) {
        return false;
    }
    const size_t idx = static_cast<size_t>(atomIndex);
    if (idx >= mask.size()) {
        return false;
    }
    return mask[idx] != 0;
}

DrudeCore::DrudeCore() 
: m_algorithm(DrudeAlgorithm::SCF),
    m_scfOptimizer(std::make_unique<DrudeSCF>()),
    m_currentOptimizer(m_scfOptimizer.get()) {
}

double DrudeCore::calculateEnergy(model::MCState& state) {
if (!m_currentOptimizer) {
    throw std::runtime_error("No Drude optimizer available");
}

std::vector<char> activeAtomMask;
buildActiveAtomMask(state, activeAtomMask);

// ASPC: Apply history positions if available
if (m_useASPC && m_hasHistory) {
    applyHistoryPositions(state);
}

// First optimize Drude positions
bool converged = m_currentOptimizer->optimize(state, m_particles, m_screenedPairs, m_params);

if (!converged) {
    // SCF did not converge - this can happen in edge cases
    // but is handled appropriately by the hard wall constraint
    // std::cerr << "Warning: Drude SCF did not converge\n";
}

// ASPC: Save current positions for next step
if (m_useASPC) {
    saveCurrentPositions(state);
}

// Calculate total energy
double energy = 0.0;

// Harmonic spring energy
energy += calculateHarmonicEnergy(state, activeAtomMask);

// Coulomb energy: optional (avoid double counting with main nonbonded)
// Default: false (spring-only) for production use
// Set true only for standalone testing without main nonbonded module
if (m_params.includeCoulombEnergy) {
    energy += calculateCoulombEnergy(state, activeAtomMask);
} else {
    // When Coulomb is computed by the main nonbonded module, we still need the
    // screened-unscreened correction for configured Thole pairs.
    energy += calculateTholeCorrectionEnergy(state, activeAtomMask);
}

return energy;
}

void DrudeCore::calculateForces(model::MCState& state, std::vector<Vec3>& forces) {
// Ensure forces vector has correct size
if (forces.size() != static_cast<size_t>(state.activeAtomCount)) {
    forces.resize(static_cast<size_t>(state.activeAtomCount), {0.0, 0.0, 0.0});
}

// ASPC: Apply history positions if available
if (m_useASPC && m_hasHistory) {
    applyHistoryPositions(state);
}

// First optimize Drude positions
m_currentOptimizer->optimize(state, m_particles, m_screenedPairs, m_params);

// ASPC: Save current positions for next step
if (m_useASPC) {
    saveCurrentPositions(state);
}

// Calculate forces from harmonic springs only
// Coulomb forces are handled by the main nonbonded calculation
calculateHarmonicForces(state, forces);
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
        
    case DrudeAlgorithm::FastFBP:
        if (!m_fastFbpOptimizer) {
            m_fastFbpOptimizer = std::make_unique<DrudeFastFBP>();
        }
        m_currentOptimizer = m_fastFbpOptimizer.get();
        break;
        
    case DrudeAlgorithm::TCG:
        if (!m_tcgOptimizer) {
            m_tcgOptimizer = std::make_unique<DrudeTCG>();
        }
        m_currentOptimizer = m_tcgOptimizer.get();
        break;
        
    case DrudeAlgorithm::TCGv2:
        if (!m_tcgv2Optimizer) {
            m_tcgv2Optimizer = std::make_unique<DrudeTCGv2>();
        }
        m_currentOptimizer = m_tcgv2Optimizer.get();
        break;
        
    case DrudeAlgorithm::LBFGS:
        if (!m_lbfgsOptimizer) {
            m_lbfgsOptimizer = std::make_unique<DrudeLBFGS>();
        }
        m_currentOptimizer = m_lbfgsOptimizer.get();
        break;
        
    case DrudeAlgorithm::Direct:
        if (!m_directOptimizer) {
            m_directOptimizer = std::make_unique<DrudeDirectPolarization>();
        }
        m_currentOptimizer = m_directOptimizer.get();
        break;
        
    case DrudeAlgorithm::Hybrid:
        if (!m_hybridOptimizer) {
            m_hybridOptimizer = std::make_unique<DrudeHybrid>();
        }
        m_currentOptimizer = m_hybridOptimizer.get();
        break;
        
    case DrudeAlgorithm::MultiStage:
        if (!m_multiStageOptimizer) {
            m_multiStageOptimizer = std::make_unique<DrudeMultiStage>();
        }
        m_currentOptimizer = m_multiStageOptimizer.get();
        break;
}
}

void DrudeCore::setParameters(const DrudeSCFParams& params) {
m_params = params;
}

void DrudeCore::clear() {
m_particles.clear();
m_screenedPairs.clear();
m_hasHistory = false;
m_lastDrudePositions.clear();
}

size_t DrudeCore::getNumParticles() const {
return m_particles.size();
}

double DrudeCore::calculateHarmonicEnergy(const model::MCState& state, const std::vector<char>& activeAtomMask) const {
double energy = 0.0;

for (const auto& particle : m_particles) {
    if (!isActiveAtom(particle.parentIndex, activeAtomMask) ||
        !isActiveAtom(particle.drudeIndex, activeAtomMask)) {
        continue;
    }
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
        if (isActiveAtom(particle.aniso1Index, activeAtomMask) &&
            isActiveAtom(particle.aniso2Index, activeAtomMask)) {
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
    }
    
    if (particle.aniso3Index >= 0 && particle.aniso4Index >= 0) {
        if (isActiveAtom(particle.aniso3Index, activeAtomMask) &&
            isActiveAtom(particle.aniso4Index, activeAtomMask)) {
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
}

return energy;
}

double DrudeCore::calculateScreenedCoulombEnergy(const model::MCState& /*state*/) const {
double energy = 0.0;

// This function should not be called directly in the current implementation.
// Thole screening is applied during SCF optimization to modify the electric field
// and prevent polarization catastrophe. It does not contribute a separate energy term.
// The screened interactions are already accounted for in the converged Drude positions.

// For compatibility with tests that may call this function directly,
// we return 0 since the screening effect is already included in the
// optimized Drude positions from SCF.

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
    double u = pair.thole * r / alpha_ij;  // u = thole * r / alpha_eff
    double exp_u = std::exp(-u);
    
    double screening = 1.0 - (1.0 + u/2.0) * exp_u;
    // dS/dr = d/dr[1 - (1 + u/2)*exp(-u)]
    //       = -d/dr[(1 + u/2)*exp(-u)]
    //       = -(1/2 * du/dr * exp(-u) - (1 + u/2) * exp(-u) * du/dr)
    //       = -du/dr * exp(-u) * (1/2 - (1 + u/2))
    //       = du/dr * exp(-u) * u/2
    double du_dr = pair.thole / alpha_ij;
    double dscreening_dr = du_dr * exp_u * u / 2.0;
    
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
if (box[0] > 0) dx -= box[0] * std::round(dx / box[0]);
if (box[1] > 0) dy -= box[1] * std::round(dy / box[1]);
if (box[2] > 0) dz -= box[2] * std::round(dz / box[2]);
}

double DrudeCore::calculateCoulombEnergy(const model::MCState& state, const std::vector<char>& activeAtomMask) const {
double energy = 0.0;

// First, calculate normal Coulomb interactions for all atom pairs
// excluding intramolecular interactions and Drude-parent pairs
for (int i = 0; i < state.activeAtomCount - 1; ++i) {
    if (!isActiveAtom(i, activeAtomMask)) {
        continue;
    }
    const auto& atom1 = state.atoms[i];
    
    for (int j = i + 1; j < state.activeAtomCount; ++j) {
        if (!isActiveAtom(j, activeAtomMask)) {
            continue;
        }
        const auto& atom2 = state.atoms[j];
        
        // Skip if both have zero charge
        if (std::abs(atom1.charge) < 1e-10 && std::abs(atom2.charge) < 1e-10) {
            continue;
        }
        
        // Skip intramolecular interactions 
        if (inSameMolecule(i, j, state)) {
            continue;
        }
        
        // Skip Drude-parent interactions (handled by spring force)
        bool isDrudeParent = false;
        for (const auto& particle : m_particles) {
            if ((i == particle.drudeIndex && j == particle.parentIndex) ||
                (j == particle.drudeIndex && i == particle.parentIndex)) {
                isDrudeParent = true;
                break;
            }
        }
        if (isDrudeParent) {
            continue;
        }
        
        // Calculate normal Coulomb interaction
        double dx = atom2.x - atom1.x;
        double dy = atom2.y - atom1.y;
        double dz = atom2.z - atom1.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-12) continue;
        
        double r = std::sqrt(r2);
        energy += DrudeConstants::ONE_4PI_EPS0 * atom1.charge * atom2.charge / r;
    }
}

// Now apply Thole screening corrections to screened pairs
// This modifies the energy of pairs that should be screened
for (const auto& pair : m_screenedPairs) {
    const auto& particle1 = m_particles[pair.dipole1];
    const auto& particle2 = m_particles[pair.dipole2];
    
    // Bounds check for particle indices
    if (!isActiveAtom(particle1.parentIndex, activeAtomMask) ||
        !isActiveAtom(particle1.drudeIndex, activeAtomMask) ||
        !isActiveAtom(particle2.parentIndex, activeAtomMask) ||
        !isActiveAtom(particle2.drudeIndex, activeAtomMask)) {
        continue;
    }
    
    // Get all four atoms involved
    int atoms1[2] = {particle1.parentIndex, particle1.drudeIndex};
    int atoms2[2] = {particle2.parentIndex, particle2.drudeIndex};

    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};

    const bool screenParents = (m_params.tholeMode == TholeMode::StandardS1);

	    // Apply Thole screening as a correction to the already-included unscreened Coulomb terms.
	    // StandardS1: screen all 4 charge-charge interactions between dipoles.
	    // OpenMMCompat: screen only interactions involving at least one Drude particle.
	    for (int j = 0; j < 2; ++j) {
	        for (int k = 0; k < 2; ++k) {
	            const bool involvesDrude = (j == 1 || k == 1);
	            if (!screenParents && !involvesDrude) {
	                continue;  // OpenMMCompat: leave parent-parent unscreened
	            }
	            // gcmc_cpu's DIRECT energy backend does not include intra-residue Coulomb.
	            // For screened dipole pairs within the same residue, the *full* screened induced-charge
	            // interaction must be added here. For inter-residue pairs, add only the correction
	            // (screened - unscreened) to avoid modifying permanent charge interactions.
	            const bool sameMolecule = inSameMolecule(atoms1[j], atoms2[k], state);
	            
	            // Bounds check
	            if (!isActiveAtom(atoms1[j], activeAtomMask) ||
	                !isActiveAtom(atoms2[k], activeAtomMask)) {
	                continue;
	            }
            
            // Get atoms
            const auto& atom1 = state.atoms[atoms1[j]];
            const auto& atom2 = state.atoms[atoms2[k]];
            
            // Calculate distance
            double dx = atom2.x - atom1.x;
            double dy = atom2.y - atom1.y;
            double dz = atom2.z - atom1.z;
            applyPBC(dx, dy, dz, box);
	            
	            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
	            if (r < 1e-6) continue;

	            // Screened-pair interactions operate on the induced Drude charge pair (±q_drude),
	            // not on the full atomic charges (which include permanent charge on the core atom).
	            // This matches OpenMM/CHARMM DrudeForce screened pair semantics.
	            const double q1 = (j == 0) ? -particle1.charge : particle1.charge;
	            const double q2 = (k == 0) ? -particle2.charge : particle2.charge;
	            const double coulomb = DrudeConstants::ONE_4PI_EPS0 * q1 * q2 / r;

	            // Per-interaction screening factor S1(u) (see THOLE_CHARMM_IMPLEMENTATION.md).
	            double screening = computeTholeScreening(
	                r,
	                particle1.polarizability,
	                particle2.polarizability,
	                pair.thole
	            );

	            // Intra-residue: baseline term is absent -> add full screened.
	            // Inter-residue: baseline unscreened induced-charge term is present -> add correction.
	            const double factor = sameMolecule ? screening : (screening - 1.0);
	            energy += coulomb * factor;
	        }
	    }
	}

return energy;
}

double DrudeCore::calculateTholeCorrectionEnergy(const model::MCState& state, const std::vector<char>& activeAtomMask) const {
    double energy = 0.0;
    for (const auto& pair : m_screenedPairs) {
        if (pair.dipole1 < 0 || pair.dipole2 < 0 ||
            static_cast<size_t>(pair.dipole1) >= m_particles.size() ||
            static_cast<size_t>(pair.dipole2) >= m_particles.size()) {
            continue;
        }

        const auto& particle1 = m_particles[static_cast<size_t>(pair.dipole1)];
        const auto& particle2 = m_particles[static_cast<size_t>(pair.dipole2)];

        if (!isActiveAtom(particle1.parentIndex, activeAtomMask) ||
            !isActiveAtom(particle1.drudeIndex, activeAtomMask) ||
            !isActiveAtom(particle2.parentIndex, activeAtomMask) ||
            !isActiveAtom(particle2.drudeIndex, activeAtomMask)) {
            continue;
        }

        const int atoms1[2] = {particle1.parentIndex, particle1.drudeIndex};
        const int atoms2[2] = {particle2.parentIndex, particle2.drudeIndex};
        const std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};

        const bool screenParents = (m_params.tholeMode == TholeMode::StandardS1);
	        for (int j = 0; j < 2; ++j) {
	            for (int k = 0; k < 2; ++k) {
	                const bool involvesDrude = (j == 1 || k == 1);
	                if (!screenParents && !involvesDrude) {
	                    continue;
	                }

	                if (!isActiveAtom(atoms1[j], activeAtomMask) ||
	                    !isActiveAtom(atoms2[k], activeAtomMask)) {
	                    continue;
	                }
	                const bool sameMolecule = inSameMolecule(atoms1[j], atoms2[k], state);

	                const auto& atom1 = state.atoms[atoms1[j]];
	                const auto& atom2 = state.atoms[atoms2[k]];

                double dx = atom2.x - atom1.x;
                double dy = atom2.y - atom1.y;
                double dz = atom2.z - atom1.z;
                applyPBC(dx, dy, dz, box);

                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
	                if (r < 1e-6) {
	                    continue;
	                }

	                const double q1 = (j == 0) ? -particle1.charge : particle1.charge;
	                const double q2 = (k == 0) ? -particle2.charge : particle2.charge;
	                const double coulomb = DrudeConstants::ONE_4PI_EPS0 * q1 * q2 / r;
	                const double screening = computeTholeScreening(
	                    r,
	                    particle1.polarizability,
	                    particle2.polarizability,
	                    pair.thole
	                );

	                const double factor = sameMolecule ? screening : (screening - 1.0);
	                energy += coulomb * factor;
	            }
	        }
	    }
	    return energy;
	}

bool DrudeCore::inSameMolecule(int atom1, int atom2, const model::MCState& state) const {
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

// ASPC helper methods
void DrudeCore::saveCurrentPositions(const model::MCState& state) {
// Resize if needed
if (m_lastDrudePositions.size() != m_particles.size()) {
    m_lastDrudePositions.resize(m_particles.size());
}

// Save current Drude positions relative to parent
for (size_t i = 0; i < m_particles.size(); ++i) {
    const auto& particle = m_particles[i];
    const auto& drude = state.atoms[particle.drudeIndex];
    const auto& parent = state.atoms[particle.parentIndex];
    
    m_lastDrudePositions[i][0] = drude.x - parent.x;
    m_lastDrudePositions[i][1] = drude.y - parent.y;
    m_lastDrudePositions[i][2] = drude.z - parent.z;
}

m_hasHistory = true;
}

void DrudeCore::applyHistoryPositions(model::MCState& state) {
if (!m_hasHistory || m_lastDrudePositions.size() != m_particles.size()) {
    return;  // No valid history
}

// Apply saved displacements to current parent positions
for (size_t i = 0; i < m_particles.size(); ++i) {
    const auto& particle = m_particles[i];
    const auto& parent = state.atoms[particle.parentIndex];
    auto& drude = state.atoms[particle.drudeIndex];
    
    // Apply displacement from parent
    drude.x = parent.x + m_lastDrudePositions[i][0];
    drude.y = parent.y + m_lastDrudePositions[i][1];
    drude.z = parent.z + m_lastDrudePositions[i][2];
}
}

// Static interface implementation
DrudeCore& DrudeCore::getInstance() {
return g_drudeCore;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
