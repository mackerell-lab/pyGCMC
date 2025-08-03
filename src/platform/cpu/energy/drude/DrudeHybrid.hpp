#pragma once

/**
 * @file DrudeHybrid.hpp
 * @brief Hybrid optimization strategy combining FastFBP with SCF
 * 
 * Uses FastFBP as a preconditioner for SCF to achieve both speed and accuracy
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include "DrudeSCF.hpp"
#include "DrudeFastFBP.hpp"
#include <memory>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Hybrid FastFBP+SCF optimizer
 * 
 * Strategy:
 * 1. Use FastFBP for initial guess (3-5 iterations)
 * 2. Switch to SCF for final convergence
 * 3. Track convergence to optimize switching point
 */
class DrudeHybrid : public DrudeOptimizer {
public:
    /**
     * @brief Hybrid strategy modes
     */
    enum class HybridMode {
        Fixed,      // Fixed FastFBP iterations before SCF
        Dynamic,    // Dynamic switching based on convergence
        Adaptive    // Adaptive with learning
    };
    
    DrudeHybrid();
    ~DrudeHybrid() = default;
    
    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;
    
    const char* getName() const override { return "Hybrid-FastFBP+SCF"; }
    
    // Configuration methods
    void setHybridMode(HybridMode mode) { hybridMode_ = mode; }
    void setFastFBPIterations(int iter) { fbpIterations_ = iter; }
    void setSwitchingThreshold(double threshold) { switchThreshold_ = threshold; }
    void setMaxSCFIterations(int iter) { maxSCFIter_ = iter; }
    
    // Get statistics
    struct HybridStats {
        int fbpIterations;      // Actual FastFBP iterations used
        int scfIterations;      // SCF iterations after switching
        double fbpTime;         // Time spent in FastFBP
        double scfTime;         // Time spent in SCF
        double switchError;     // Error when switching to SCF
        bool converged;         // Final convergence status
    };
    
    const HybridStats& getStats() const { return stats_; }
    
private:
    // Sub-optimizers
    std::unique_ptr<DrudeFastFBP> fastFBP_;
    std::unique_ptr<DrudeSCF> scf_;
    
    // Configuration
    HybridMode hybridMode_ = HybridMode::Dynamic;
    int fbpIterations_ = 5;          // Default FastFBP iterations
    double switchThreshold_ = 0.01;   // Switch when error < threshold
    int maxSCFIter_ = 20;            // Max SCF iterations after switch
    
    // Statistics
    mutable HybridStats stats_;
    
    // Helper methods
    double calculateError(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles
    );
    
    bool shouldSwitchToSCF(
        double currentError,
        int fbpIterations,
        double convergenceRate
    );
    
    void configureFastFBP(const DrudeSCFParams& params);
    void configureSCF(const DrudeSCFParams& params);
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc