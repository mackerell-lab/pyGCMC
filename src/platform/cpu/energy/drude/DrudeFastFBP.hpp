#pragma once

/**
 * @file DrudeFastFBP.hpp
 * @brief Fast Force Balance Predictor for Drude oscillators
 *
 * This implementation targets 5% accuracy for GCMC applications
 * by using direct force balance with minimal iterations
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Fast FBP implementation for GCMC
 *
 * Algorithm:
 * 1. Direct FBP: r_drude = r_parent - (q/k) * E_fixed
 * 2. Iterative correction for Drude-Drude interactions (1-3 iterations)
 * 3. Distance-based cutoff for efficiency
 *
 * Target: 5% accuracy, 10-30x speedup over SCF
 */
class DrudeFastFBP : public DrudeOptimizer {
public:
    /**
     * @brief Iteration mode for FastFBP
     */
    enum class IterationMode {
        Fixed,      // Use fixed number of iterations
        Dynamic,    // Dynamic iterations based on convergence
        Adaptive    // Adaptive with parameter adjustment
    };

    DrudeFastFBP() = default;
    ~DrudeFastFBP() = default;

    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;

    const char* getName() const override { return "FastFBP"; }

    /**
     * @brief Set number of FBP iterations (default: 2)
     */
    void setIterations(int iter) { fbpIterations_ = iter; }

    /**
     * @brief Set interaction cutoff for Drude-Drude (default: 1.0 nm)
     */
    void setCutoff(double cutoff) { drudeCutoff_ = cutoff; }

    /**
     * @brief Enable/disable Drude-Drude interactions (default: true)
     */
    void setIncludeDrudeDrude(bool include) { includeDrudeDrude_ = include; }

    /**
     * @brief Set iteration mode
     */
    void setIterationMode(IterationMode mode) { iterMode_ = mode; }

    /**
     * @brief Set convergence tolerance for dynamic mode
     */
    void setConvergenceTolerance(double tol) { convTolerance_ = tol; }

    /**
     * @brief Set maximum iterations for dynamic mode
     */
    void setMaxIterations(int max) { maxIterations_ = max; }

    /**
     * @brief Enable adaptive parameter adjustment
     */
    void setAdaptiveMode(bool enable) { adaptiveMode_ = enable; }

    /**
     * @brief Get convergence statistics
     */
    struct ConvergenceStats {
        int actualIterations;
        double finalError;
        double convergenceRate;
        bool converged;
    };

    const ConvergenceStats& getStats() const { return stats_; }

private:
    // Fixed mode parameters
    int fbpIterations_ = 3;           // Number of FBP iterations (optimized for speed)
    double drudeCutoff_ = 0.8;        // Cutoff for Drude-Drude interactions (nm)
    bool includeDrudeDrude_ = true;   // Include Drude-Drude interactions
    double dampingFactor_ = 0.7;      // Damping for stability

    // Dynamic mode parameters
    IterationMode iterMode_ = IterationMode::Fixed;
    double convTolerance_ = 1e-5;     // Convergence tolerance (nm)
    int maxIterations_ = 20;          // Maximum iterations
    int minIterations_ = 2;           // Minimum iterations
    bool adaptiveMode_ = false;       // Enable adaptive parameter adjustment

    // Convergence statistics
    mutable ConvergenceStats stats_;

    // Adaptive parameters
    struct AdaptiveParams {
        double dampingFactor = 0.7;
        double cutoff = 0.8;
        double convergenceRate = 0.0;
        int stagnationCount = 0;
    };
    AdaptiveParams adaptiveParams_;

    /**
     * @brief Compute electric field at Drude position from fixed charges
     */
    Vec3 computeFixedField(
        const model::MCState& state,
        const Vec3& drudePos,
        int drudeIdx,
        int parentIdx
    );

    /**
     * @brief Compute electric field from other Drude particles
     */
    Vec3 computeDrudeField(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const Vec3& drudePos,
        int currentDrudeIdx
    );

    /**
     * @brief Apply distance constraints if needed
     */
    void applyConstraints(
        Vec3& drudePos,
        const Vec3& parentPos,
        double maxDistance
    );

    /**
     * @brief Check if two atoms are in the same molecule
     */
    bool inSameMolecule(
        int atom1,
        int atom2,
        const model::MCState& state
    );

    /**
     * @brief Check convergence for dynamic iteration mode
     */
    bool checkConvergence(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<Vec3>& previousPositions,
        double dispTolerance,
        double forceTolerance
    );

    /**
     * @brief Save current Drude positions
     */
    void saveDrudePositions(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        std::vector<Vec3>& positions
    );

    /**
     * @brief Calculate net force on Drude particle
     */
    Vec3 calculateNetForce(
        const model::MCState& state,
        const DrudeParticle& particle,
        size_t particleIndex,
        const std::vector<DrudeParticle>& allParticles
    );

    /**
     * @brief Configure parameters based on system size
     */
    void configureForSystem(
        const model::MCState& state,
        size_t nParticles
    );

    /**
     * @brief Adjust adaptive parameters based on convergence rate
     */
    void adjustAdaptiveParameters(
        double convergenceRate
    );

    /**
     * @brief Calculate convergence rate from history
     */
    double calculateConvergenceRate(
        const std::vector<double>& errorHistory
    );
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
