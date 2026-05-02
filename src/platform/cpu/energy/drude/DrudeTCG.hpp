#pragma once

/**
 * @file DrudeTCG.hpp
 * @brief Truncated Conjugate Gradient optimizer for Drude oscillators
 *
 * Based on:
 * - Aviat et al. (2017) "Truncated Conjugate Gradient: An Optimal Strategy
 *   for the Analytical Evaluation of the Many-Body Polarization Energy and
 *   Forces in Molecular Simulations"
 * - Simmonett et al. (2021) "Efficient treatment of induced dipoles"
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include <vector>
#include <unordered_map>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Truncated Conjugate Gradient optimizer
 *
 * TCG performs a fixed number of CG iterations (typically 3-5) with
 * analytical force computation. This gives:
 * - Guaranteed energy conservation in MD
 * - 2-3x speedup over full SCF
 * - ~1-5% error in dipole moments
 *
 * Algorithm:
 * 1. Initialize dipoles to zero (or from predictor)
 * 2. Perform N_tcg CG iterations
 * 3. Use Chebyshev polynomial extrapolation for final dipole
 * 4. Compute analytical forces with chain rule
 */
class DrudeTCG : public DrudeOptimizer {
public:
    DrudeTCG() = default;
    ~DrudeTCG() = default;

    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;

    const char* getName() const override { return "TCG"; }

    /**
     * @brief Set number of CG iterations (default: 4)
     */
    void setIterations(int iter) { tcgIterations_ = iter; }

    /**
     * @brief Enable Chebyshev extrapolation (default: true)
     */
    void setUseChebyshev(bool use) { useChebyshev_ = use; }

    /**
     * @brief Set Chebyshev order (default: 3)
     */
    void setChebyshevOrder(int order) { chebyshevOrder_ = order; }

private:
    int tcgIterations_ = 4;        // Number of CG iterations
    bool useChebyshev_ = true;     // Use Chebyshev extrapolation
    int chebyshevOrder_ = 3;       // Chebyshev polynomial order

    // Working arrays (allocated once, reused)
    std::vector<Vec3> residuals_;   // Current residuals (r)
    std::vector<Vec3> directions_;  // Search directions (p)
    std::vector<Vec3> Ap_;         // Matrix-vector product
    std::vector<Vec3> oldResiduals_; // Previous residuals

    // Screening map for current optimization
    std::unordered_map<uint64_t, double> screeningMap_;

    // Chebyshev extrapolation history
    std::vector<std::vector<Vec3>> dipoleHistory_;

    /**
     * @brief Compute electric field at parent positions
     * @param state Molecular state
     * @param particles Drude particles
     * @param fields Output field vectors
     */
    void computeFieldsAtParents(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        std::vector<Vec3>& fields
    );

    /**
     * @brief Update Drude positions based on displacements
     * @param state Molecular state (will be modified)
     * @param particles Drude particles
     * @param displacements Displacement vectors from parent
     */
    void updateDrudePositions(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<Vec3>& displacements
    );

    /**
     * @brief Compute induced fields from displacements
     * @param state Molecular state
     * @param particles Drude particles
     * @param displacements Current displacements
     * @param fields Output induced fields at Drude positions
     */
    void computeInducedFieldsFromDisplacements(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<Vec3>& displacements,
        std::vector<Vec3>& fields
    );

    /**
     * @brief Compute fields from all charges
     * @param state Molecular state
     * @param particles Drude particles
     * @param fields Output field vectors
     */
    void computeFields(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        std::vector<Vec3>& fields
    );

    /**
     * @brief Check if two atoms are excluded
     */
    bool isExcluded(int atom1, int atom2, const model::MCState& state);
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
