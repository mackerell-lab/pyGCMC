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
    
private:
    int fbpIterations_ = 3;           // Number of FBP iterations (optimized for speed)
    double drudeCutoff_ = 0.8;        // Cutoff for Drude-Drude interactions (nm)
    bool includeDrudeDrude_ = true;   // Include Drude-Drude interactions
    double dampingFactor_ = 0.7;      // Damping for stability
    
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
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc