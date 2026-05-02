#pragma once

/**
 * @file DrudeSCF.hpp
 * @brief Self-Consistent Field (SCF) optimizer for Drude oscillators
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief SCF optimizer for Drude positions
 *
 * Iteratively optimizes Drude positions to minimize total energy
 * using self-consistent field approach with adaptive damping.
 */
class DrudeSCF : public DrudeOptimizer {
public:
    DrudeSCF() = default;
    ~DrudeSCF() = default;

    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;

    const char* getName() const override { return "SCF"; }

    /**
     * @brief Get the number of iterations from last optimization
     * @return Number of iterations used
     */
    int getIterationCount() const { return m_lastIterationCount; }

private:
    mutable int m_lastIterationCount = 0;

    static void buildActiveAtomMask(const model::MCState& state, std::vector<char>& mask);
    static void buildActiveAtomIndices(const std::vector<char>& mask, std::vector<int>& indices);
    static bool isActiveAtom(int atomIndex, const std::vector<char>& mask);

    /**
     * @brief Calculate electric field at Drude particles
     *
     * Includes contributions from:
     * 1. External charges (atoms)
     * 2. Other Drude particles (with Thole screening)
     */
    void calculateElectricField(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        std::vector<Vec3>& electricField,
        const std::vector<int>& activeAtoms,
        const std::vector<char>& activeAtomMask
    ) const;

    /**
     * @brief Calculate field from external charges
     */
    void calculateExternalField(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        std::vector<Vec3>& electricField,
        const std::vector<int>& activeAtoms
    ) const;

    /**
     * @brief Calculate field from other Drude particles
     */
    void calculateInducedField(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        std::vector<Vec3>& electricField,
        const std::vector<char>& activeAtomMask
    ) const;

    /**
     * @brief Update Drude positions based on electric field
     * @return Maximum displacement
     */
    double updateDrudePositions(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<Vec3>& electricField,
        double dampingFactor,
        double maxDrudeDistance,
        const std::vector<char>& activeAtomMask
    ) const;

    /**
     * @brief Calculate force on each Drude particle
     */
    void calculateDrudeForces(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<Vec3>& electricField,
        std::vector<Vec3>& forces,
        const std::vector<char>& activeAtomMask
    ) const;

    /**
     * @brief Apply periodic boundary conditions
     */
    void applyPBC(double& dx, double& dy, double& dz, const std::array<double, 3>& box) const;

    /**
     * @brief Check if two atoms are in the same molecule
     */
    bool inSameMolecule(int atom1, int atom2, const model::MCState& state) const;
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
