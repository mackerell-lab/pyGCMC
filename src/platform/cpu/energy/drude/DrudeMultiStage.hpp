#pragma once

/**
 * @file DrudeMultiStage.hpp
 * @brief Multi-stage optimization strategy: Direct → FastFBP → TCG → SCF
 *
 * This implements the optimal hybrid strategy identified through analysis
 * Enhanced based on literature review showing TCG-3 achieves <1% error with ~15x speedup
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include "DrudeDirectPolarization.hpp"
#include "DrudeFastFBP.hpp"
#include "DrudeTCG.hpp"
#include "DrudeSCF.hpp"
#include <memory>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Multi-stage optimizer combining Direct, FastFBP, TCG and SCF
 *
 * Optimization stages:
 * 1. Direct polarization - instant physical guess (0 iterations)
 * 2. FastFBP - rapid improvement (3-10 iterations)
 * 3. TCG - efficient refinement (3-5 iterations, <1% error)
 * 4. SCF - fine convergence (adaptive iterations, optional)
 *
 * Based on literature:
 * - TCG-3 achieves <1% error with ~15x speedup (Aviat et al. 2017)
 * - Direct polarization effective for GCMC (Drew & Gilson 2025)
 */
class DrudeMultiStage : public DrudeOptimizer {
public:
    /**
     * @brief Configuration for multi-stage optimization
     */
    struct MultiStageConfig {
        // FastFBP stage
        int minFBPIterations = 3;
        int maxFBPIterations = 10;
        double fbpCutoffFactor = 0.8;  // Fraction of full cutoff

        // TCG stage
        int tcgIterations = 3;           // TCG-3 recommended by literature
        bool enableTCG = true;           // Enable TCG stage
        double tcgErrorThreshold = 0.01; // Switch to TCG when error < this

        // Switching criteria
        double switchToSCFError = 0.001;  // Switch when error < this (tighter)
        double switchToSCFRate = 0.05;    // Switch when convergence rate < this

        // SCF stage
        double scfTolerance = 0.001;      // Final tolerance (tighter)
        int maxSCFIterations = 20;        // Fewer needed after TCG
        bool requireSCF = false;          // SCF optional after TCG

        // Adaptive parameters
        bool adaptiveMode = true;
        double densityThreshold = 1.0;   // g/cm³
        double polarizabilityThreshold = 0.001;  // nm³
    };

    /**
     * @brief Statistics for performance analysis
     */
    struct MultiStageStats {
        // Timing
        double directTime = 0.0;
        double fbpTime = 0.0;
        double tcgTime = 0.0;
        double scfTime = 0.0;

        // Iterations
        int fbpIterations = 0;
        int tcgIterations = 0;
        int scfIterations = 0;

        // Errors at each stage
        double errorAfterDirect = 0.0;
        double errorAfterFBP = 0.0;
        double errorAfterTCG = 0.0;
        double finalError = 0.0;

        // System characteristics
        double systemDensity = 0.0;
        double avgPolarizability = 0.0;

        bool converged = false;
    };

    DrudeMultiStage();
    ~DrudeMultiStage() = default;

    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;

    const char* getName() const override { return "MultiStage-Direct+FastFBP+TCG+SCF"; }

    // Configuration
    void setConfig(const MultiStageConfig& config) { config_ = config; }
    const MultiStageConfig& getConfig() const { return config_; }

    // Statistics
    const MultiStageStats& getStats() const { return stats_; }

    // Enable/disable stages
    void enableDirectStage(bool enable) { useDirectStage_ = enable; }
    void enableFBPStage(bool enable) { useFBPStage_ = enable; }
    void enableTCGStage(bool enable) { useTCGStage_ = enable; }
    void enableSCFStage(bool enable) { useSCFStage_ = enable; }

private:
    // Sub-optimizers
    std::unique_ptr<DrudeDirectPolarization> direct_;
    std::unique_ptr<DrudeFastFBP> fastFBP_;
    std::unique_ptr<DrudeTCG> tcg_;
    std::unique_ptr<DrudeSCF> scf_;

    // Configuration
    MultiStageConfig config_;
    mutable MultiStageStats stats_;

    // Stage control
    bool useDirectStage_ = true;
    bool useFBPStage_ = true;
    bool useTCGStage_ = true;
    bool useSCFStage_ = true;

    // Helper methods
    void analyzeSystem(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles
    );

    void adaptConfiguration(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles
    );

    double calculateError(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles
    ) const;

    int determineFBPIterations(
        double density,
        double polarizability,
        size_t systemSize
    ) const;

    bool shouldSwitchToSCF(
        double currentError,
        double convergenceRate,
        int fbpIterations
    ) const;
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
