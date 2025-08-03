#pragma once

/**
 * @file DrudeTCGv2.hpp
 * @brief Improved Truncated Conjugate Gradient optimizer for Drude oscillators
 * 
 * Improvements:
 * - Correct Drude-Drude interaction handling
 * - Diagonal preconditioning
 * - Better convergence monitoring
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include <vector>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Improved TCG implementation with better accuracy
 */
class DrudeTCGv2 : public DrudeOptimizer {
public:
    DrudeTCGv2() = default;
    ~DrudeTCGv2() = default;
    
    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;
    
    const char* getName() const override { return "TCGv2"; }
    
    /**
     * @brief Set number of CG iterations (default: 5)
     */
    void setIterations(int iter) { tcgIterations_ = iter; }
    
    /**
     * @brief Enable diagonal preconditioning (default: true)
     */
    void setUsePreconditioning(bool use) { usePreconditioning_ = use; }
    
    /**
     * @brief Enable debug output (default: false)
     */
    void setDebugMode(bool debug) { debugMode_ = debug; }
    
private:
    int tcgIterations_ = 5;           // Number of CG iterations
    bool usePreconditioning_ = true;  // Use diagonal preconditioner
    bool debugMode_ = false;          // Print debug info
    
    // Working arrays
    std::vector<Vec3> residuals_;     // Current residuals (r)
    std::vector<Vec3> directions_;    // Search directions (p)
    std::vector<Vec3> Ap_;           // Matrix-vector product
    std::vector<Vec3> z_;            // Preconditioned residual
    std::vector<double> precond_;    // Diagonal preconditioner
    
    /**
     * @brief Compute electric field at parent positions from fixed charges
     */
    void computeFieldsAtParents(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        std::vector<Vec3>& fields
    );
    
    /**
     * @brief Apply matrix A = I - (q/k)*T to displacement vector
     * @param state Current state with updated Drude positions
     * @param particles Drude particles
     * @param displacements Input displacement vector
     * @param result Output: A*displacements
     */
    void applyMatrix(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<Vec3>& displacements,
        std::vector<Vec3>& result
    );
    
    /**
     * @brief Compute induced fields at current Drude positions
     */
    void computeInducedFields(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        std::vector<Vec3>& fields
    );
    
    /**
     * @brief Update Drude positions based on displacements
     */
    void updateDrudePositions(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<Vec3>& displacements
    );
    
    /**
     * @brief Build diagonal preconditioner
     */
    void buildPreconditioner(
        const std::vector<DrudeParticle>& particles
    );
    
    /**
     * @brief Apply preconditioner: z = M^{-1} * r
     */
    void applyPreconditioner(
        const std::vector<Vec3>& residuals,
        std::vector<Vec3>& z
    );
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc