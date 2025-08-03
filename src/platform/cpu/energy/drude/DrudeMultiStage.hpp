#pragma once

/**
 * @file DrudeMultiStage.hpp
 * @brief Multi-stage optimization strategy: Direct → FastFBP → SCF
 * 
 * This implements the optimal hybrid strategy identified through analysis
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include "DrudeDirectPolarization.hpp"
#include "DrudeFastFBP.hpp"
#include "DrudeSCF.hpp"
#include <memory>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Multi-stage optimizer combining Direct, FastFBP, and SCF
 * 
 * Optimization stages:
 * 1. Direct polarization - instant physical guess (0 iterations)
 * 2. FastFBP - rapid improvement (3-10 iterations)
 * 3. SCF - fine convergence (adaptive iterations)
 * 
 * This combination provides optimal speed-accuracy tradeoff
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
        
        // Switching criteria
        double switchToSCFError = 0.01;  // Switch when error < this
        double switchToSCFRate = 0.1;    // Switch when convergence rate < this
        
        // SCF stage
        double scfTolerance = 0.01;      // Final tolerance
        int maxSCFIterations = 50;
        
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
        double scfTime = 0.0;
        
        // Iterations
        int fbpIterations = 0;
        int scfIterations = 0;
        
        // Errors at each stage
        double errorAfterDirect = 0.0;
        double errorAfterFBP = 0.0;
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
    
    const char* getName() const override { return "MultiStage-Direct+FastFBP+SCF"; }
    
    // Configuration
    void setConfig(const MultiStageConfig& config) { config_ = config; }
    const MultiStageConfig& getConfig() const { return config_; }
    
    // Statistics
    const MultiStageStats& getStats() const { return stats_; }
    
    // Enable/disable stages
    void enableDirectStage(bool enable) { useDirectStage_ = enable; }
    void enableFBPStage(bool enable) { useFBPStage_ = enable; }
    
private:
    // Sub-optimizers
    std::unique_ptr<DrudeDirectPolarization> direct_;
    std::unique_ptr<DrudeFastFBP> fastFBP_;
    std::unique_ptr<DrudeSCF> scf_;
    
    // Configuration
    MultiStageConfig config_;
    mutable MultiStageStats stats_;
    
    // Stage control
    bool useDirectStage_ = true;
    bool useFBPStage_ = true;
    
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