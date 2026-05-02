/**
 * @file DrudeDirectPolarization.cpp
 * @brief Implementation of direct polarization approximation
 */

#include "DrudeDirectPolarization.hpp"
#include "../common/EnergyConstants.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

bool DrudeDirectPolarization::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    // Unused parameters
    (void)screenedPairs;
    (void)params;

    if (particles.empty()) return true;

    // Step 1: Compute permanent fields at parent positions
    std::vector<Vec3> permanentFields(particles.size());
    computePermanentFields(state, particles, permanentFields);

    // Step 2: Apply direct polarization
    // For Drude model: d = (q/k) * E_permanent
    // This ignores induced-induced interactions completely
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parent = state.atoms[particle.parentIndex];
        auto& drude = state.atoms[particle.drudeIndex];

        // Direct displacement from permanent field
        double factor = particle.charge / particle.kSpring;

        // New Drude position
        drude.x = parent.x + factor * permanentFields[i][0];
        drude.y = parent.y + factor * permanentFields[i][1];
        drude.z = parent.z + factor * permanentFields[i][2];

        // Apply constraints if needed
        if (params.enableHardWall && params.maxDrudeDistance > 0) {
            double dx = drude.x - parent.x;
            double dy = drude.y - parent.y;
            double dz = drude.z - parent.z;
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > params.maxDrudeDistance * params.maxDrudeDistance) {
                double r = std::sqrt(r2);
                double scale = params.maxDrudeDistance / r;
                drude.x = parent.x + scale * dx;
                drude.y = parent.y + scale * dy;
                drude.z = parent.z + scale * dz;
            }
        }
    }

    return true;  // Always "converges" in one step
}

void DrudeDirectPolarization::computePermanentFields(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    std::vector<Vec3>& fields
) {
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const auto& box = state.info.box;
    const double halfBox[3] = {box[0] * 0.5, box[1] * 0.5, box[2] * 0.5};

    // Initialize fields
    for (auto& field : fields) {
        field[0] = field[1] = field[2] = 0.0;
    }

    // Compute field at each parent position from permanent charges only
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parent = state.atoms[particle.parentIndex];

        // Sum over all atoms
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self
            if (j == particle.parentIndex) continue;

            // Skip all Drude particles (we want permanent charges only)
            bool isDrude = false;
            for (const auto& p : particles) {
                if (j == p.drudeIndex) {
                    isDrude = true;
                    break;
                }
            }
            if (isDrude) continue;

            // Skip intramolecular interactions
            if (inSameMolecule(particle.parentIndex, j, state)) continue;

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

bool DrudeDirectPolarization::inSameMolecule(
    int atom1,
    int atom2,
    const model::MCState& state
) {
    // Check if we have residues defined
    if (state.activeResidueCount <= 0 || state.residues.empty()) {
        return false;  // No residues defined
    }

    // Find which residue each atom belongs to
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (i >= static_cast<int>(state.residues.size())) {
            break;
        }

        const auto& res = state.residues[i];
        bool atom1InRes = (atom1 >= res.atomStart &&
                          atom1 < res.atomStart + res.atomCount);
        bool atom2InRes = (atom2 >= res.atomStart &&
                          atom2 < res.atomStart + res.atomCount);

        if (atom1InRes && atom2InRes) {
            return true;  // Both atoms in same residue
        }
    }

    return false;  // Atoms in different residues
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
