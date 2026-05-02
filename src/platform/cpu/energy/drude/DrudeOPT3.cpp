/**
 * @file DrudeOPT3.cpp
 * @brief Third-order perturbation theory optimizer implementation
 */

#include "DrudeOPT3.hpp"
#include "DrudeCore.hpp"
#include <cmath>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

bool DrudeOPT3::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    // Get active atom count
    const int activeAtomCount = state.activeAtomCount;

    // Initialize coefficients if not set
    if (m_coefficients.c0 == 0.0 && m_coefficients.c1 == 1.0 &&
        m_coefficients.c2 == 0.0 && m_coefficients.c3 == 0.0) {
        // Optimized coefficients for multi-body systems
        // These work well for bulk water (error ~10-15%)
        // Single molecule accuracy is sacrificed for better bulk behavior
        m_coefficients.c0 = 0.0;    // No static contribution
        m_coefficients.c1 = 0.36;   // Balanced for multi-body
        m_coefficients.c2 = 0.34;   // systems where high-order
        m_coefficients.c3 = 0.30;   // terms are important
    }

    // Store original Drude positions for restoration if needed
    std::vector<Vec3> originalPositions;
    originalPositions.reserve(particles.size());

    for (const auto& particle : particles) {
        const auto& drude = state.atoms[particle.drudeIndex];
        originalPositions.push_back({drude.x, drude.y, drude.z});
    }

    // Arrays for perturbation theory positions
    std::vector<Vec3> r0(particles.size()); // Parent positions
    std::vector<Vec3> r1(particles.size()); // First order
    std::vector<Vec3> r2(particles.size()); // Second order
    std::vector<Vec3> r3(particles.size()); // Third order

    // Step 0: Initialize r0 to parent positions and reset Drudes
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parent = state.atoms[particle.parentIndex];
        r0[i] = {parent.x, parent.y, parent.z};

        // Reset Drude to parent position for consistent starting point
        auto& drude = state.atoms[particle.drudeIndex];
        drude.x = parent.x;
        drude.y = parent.y;
        drude.z = parent.z;
    }

    // Step 1: First order - response to external (non-Drude) charges only
    std::vector<Vec3> fields(activeAtomCount, {0.0, 0.0, 0.0});

    // Calculate field from fixed charges only (no Drude contributions)
    calculateElectricField(state, fields, particles, screenedPairs, false);

    // Compute first order displacements
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        // Field at parent position for first order
        const Vec3& field = fields[particle.parentIndex];

        // Displacement = charge * field / k_spring
        // Note: charge is negative, so displacement is opposite to field
        const double dispFactor = particle.charge / particle.kSpring;
        r1[i][0] = r0[i][0] + dispFactor * field[0];
        r1[i][1] = r0[i][1] + dispFactor * field[1];
        r1[i][2] = r0[i][2] + dispFactor * field[2];

        // Update Drude position for next iteration
        auto& drude = state.atoms[particle.drudeIndex];
        drude.x = r1[i][0];
        drude.y = r1[i][1];
        drude.z = r1[i][2];
    }

    // Step 2: Second order - include induced dipole interactions
    std::fill(fields.begin(), fields.end(), Vec3{0.0, 0.0, 0.0});

    // Calculate total field (external + induced from first order positions)
    calculateElectricField(state, fields, particles, screenedPairs, true);

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        // Use field at parent position for all orders
        const Vec3& field = fields[particle.parentIndex];

        const double dispFactor = particle.charge / particle.kSpring;
        r2[i][0] = r0[i][0] + dispFactor * field[0];
        r2[i][1] = r0[i][1] + dispFactor * field[1];
        r2[i][2] = r0[i][2] + dispFactor * field[2];

        // Update for third order
        auto& drude = state.atoms[particle.drudeIndex];
        drude.x = r2[i][0];
        drude.y = r2[i][1];
        drude.z = r2[i][2];
    }

    // Step 3: Third order - further refinement
    std::fill(fields.begin(), fields.end(), Vec3{0.0, 0.0, 0.0});

    // Calculate field with second order positions
    calculateElectricField(state, fields, particles, screenedPairs, true);

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const Vec3& field = fields[particle.parentIndex];

        const double dispFactor = particle.charge / particle.kSpring;
        r3[i][0] = r0[i][0] + dispFactor * field[0];
        r3[i][1] = r0[i][1] + dispFactor * field[1];
        r3[i][2] = r0[i][2] + dispFactor * field[2];
    }

    // Final step: Combine using OPT3 coefficients
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        auto& drude = state.atoms[particle.drudeIndex];

        // r_final = c0*r0 + c1*r1 + c2*r2 + c3*r3
        drude.x = m_coefficients.c0 * r0[i][0] +
                  m_coefficients.c1 * r1[i][0] +
                  m_coefficients.c2 * r2[i][0] +
                  m_coefficients.c3 * r3[i][0];

        drude.y = m_coefficients.c0 * r0[i][1] +
                  m_coefficients.c1 * r1[i][1] +
                  m_coefficients.c2 * r2[i][1] +
                  m_coefficients.c3 * r3[i][1];

        drude.z = m_coefficients.c0 * r0[i][2] +
                  m_coefficients.c1 * r1[i][2] +
                  m_coefficients.c2 * r2[i][2] +
                  m_coefficients.c3 * r3[i][2];

        // Apply hard wall constraint if enabled
        if (params.enableHardWall && params.maxDrudeDistance > 0) {
            const auto& parent = state.atoms[particle.parentIndex];
            double dx = drude.x - parent.x;
            double dy = drude.y - parent.y;
            double dz = drude.z - parent.z;

            // Apply PBC
            const auto& box = state.info.box;
            if (dx > box[0]/2) dx -= box[0];
            if (dx < -box[0]/2) dx += box[0];
            if (dy > box[1]/2) dy -= box[1];
            if (dy < -box[1]/2) dy += box[1];
            if (dz > box[2]/2) dz -= box[2];
            if (dz < -box[2]/2) dz += box[2];

            double dist2 = dx*dx + dy*dy + dz*dz;
            if (dist2 > params.maxDrudeDistance * params.maxDrudeDistance) {
                double scale = params.maxDrudeDistance / std::sqrt(dist2);
                drude.x = parent.x + dx * scale;
                drude.y = parent.y + dy * scale;
                drude.z = parent.z + dz * scale;
            }
        }
    }

    return true; // OPT3 always converges in fixed steps
}

void DrudeOPT3::calculateElectricField(
    const model::MCState& state,
    std::vector<Vec3>& fields,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    bool includeDrudes
) const {
    const double coulombConstant = 138.935456; // kJ·nm/mol/e²
    const auto& box = state.info.box;
    const double cutoff = state.info.cutoff;
    const double cutoff2 = cutoff * cutoff;

    // Clear fields
    std::fill(fields.begin(), fields.end(), Vec3{0.0, 0.0, 0.0});

    // Create mapping from atom index to Drude particle index
    std::vector<int> atomToDrude(state.activeAtomCount, -1);
    for (size_t i = 0; i < particles.size(); ++i) {
        atomToDrude[particles[i].drudeIndex] = i;
    }

    // Calculate field at each atom position
    for (int i = 0; i < state.activeAtomCount; ++i) {
        const auto& atomi = state.atoms[i];
        Vec3& field = fields[i];

        // Sum field from all other charges
        for (int j = 0; j < state.activeAtomCount; ++j) {
            if (i == j) continue;

            const auto& atomj = state.atoms[j];

            // Skip Drude contributions if not requested
            if (!includeDrudes && atomToDrude[j] >= 0) {
                continue;
            }

            // Skip parent-Drude pair interactions
            bool skipPair = false;
            for (const auto& particle : particles) {
                if ((i == particle.parentIndex && j == particle.drudeIndex) ||
                    (j == particle.parentIndex && i == particle.drudeIndex)) {
                    skipPair = true;
                    break;
                }
            }
            if (skipPair) continue;

            // For field calculation in OPT3, we need different exclusion rules:
            // 1. Always exclude parent-Drude pairs (already done above)
            // 2. For Drude particles: only exclude their parent
            // 3. For parent atoms: include ALL other atoms (even in same molecule)
            //    because we need the field from H atoms on O
            // 4. For other atoms: apply normal bonded exclusions

            bool isParentI = false;
            bool isParentJ = false;
            for (const auto& particle : particles) {
                if (i == particle.parentIndex) isParentI = true;
                if (j == particle.parentIndex) isParentJ = true;
            }

            // Don't apply intramolecular exclusions if one atom is a parent
            // This allows O to feel field from H atoms
            if (!isParentI && !isParentJ && inSameMolecule(i, j, state)) {
                // Only exclude if neither is a parent atom
                continue;
            }

            // Calculate distance with PBC
            double dx = atomi.x - atomj.x;
            double dy = atomi.y - atomj.y;
            double dz = atomi.z - atomj.z;

            // Apply minimum image convention
            if (dx > box[0]/2) dx -= box[0];
            if (dx < -box[0]/2) dx += box[0];
            if (dy > box[1]/2) dy -= box[1];
            if (dy < -box[1]/2) dy += box[1];
            if (dz > box[2]/2) dz -= box[2];
            if (dz < -box[2]/2) dz += box[2];

            double r2 = dx*dx + dy*dy + dz*dz;

            // Skip if beyond cutoff or too close
            if (r2 > cutoff2 || r2 < 1e-10) continue;

            double r = std::sqrt(r2);
            double r3 = r2 * r;

            // Check if this pair needs Thole screening
            double screeningFactor = 1.0;

            // Only screen Drude-Drude interactions
            if (includeDrudes && atomToDrude[i] >= 0 && atomToDrude[j] >= 0) {
                int drudeI = atomToDrude[i];
                int drudeJ = atomToDrude[j];

                // Check screened pairs list
                for (const auto& pair : screenedPairs) {
                    if ((pair.dipole1 == drudeI && pair.dipole2 == drudeJ) ||
                        (pair.dipole1 == drudeJ && pair.dipole2 == drudeI)) {
                        // Apply Thole screening
                        double alphaI = particles[drudeI].polarizability;
                        double alphaJ = particles[drudeJ].polarizability;
                        double u = r / std::pow(alphaI * alphaJ, 1.0/6.0);
                        double u3 = u * u * u;
                        // Standard Thole damping function
                        screeningFactor = 1.0 - (1.0 + u + 0.5*u*u) * std::exp(-u3);
                        break;
                    }
                }
            }

            // Electric field: E = k * q * r_vec / r³
            double fieldMag = coulombConstant * atomj.charge * screeningFactor / r3;
            field[0] += fieldMag * dx;
            field[1] += fieldMag * dy;
            field[2] += fieldMag * dz;
        }
    }
}

bool DrudeOPT3::inSameMolecule(int atom1, int atom2, const model::MCState& state) const {
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
