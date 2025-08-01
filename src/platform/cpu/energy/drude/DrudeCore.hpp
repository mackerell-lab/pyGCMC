#pragma once

/**
 * @file DrudeCore.hpp
 * @brief Core implementation of Drude oscillator calculations
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include "DrudeSCF.hpp"
#include "DrudeOPT3.hpp"
#include "DrudeFBP.hpp"
#include "DrudeLBFGS.hpp"
#include <memory>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Core implementation of Drude force calculations
 * 
 * This class manages Drude particles, implements energy/force calculations,
 * and delegates optimization to specific algorithm implementations.
 */
class DrudeCore : public DrudeInterface {
public:
    DrudeCore();
    ~DrudeCore() = default;
    
    // DrudeInterface implementation
    double calculateEnergy(model::MCState& state) override;
    void calculateForces(model::MCState& state, std::vector<Vec3>& forces) override;
    int addParticle(const DrudeParticle& particle) override;
    void addScreenedPair(const ScreenedPair& pair) override;
    void setAlgorithm(DrudeAlgorithm algorithm) override;
    void setParameters(const DrudeSCFParams& params) override;
    void clear() override;
    size_t getNumParticles() const override;
    
    // Static instance access for global interface
    static DrudeCore& getInstance();
    
private:
    // Energy calculation components
    double calculateHarmonicEnergy(const model::MCState& state) const;
    double calculateScreenedCoulombEnergy(const model::MCState& state) const;
    double calculateCoulombEnergy(const model::MCState& state) const;
    
    // Force calculation components
    void calculateHarmonicForces(const model::MCState& state, 
                                std::vector<Vec3>& forces) const;
    void calculateScreenedCoulombForces(const model::MCState& state,
                                       std::vector<Vec3>& forces) const;
    
    // Utility functions
    void applyPBC(double& dx, double& dy, double& dz, const std::array<double, 3>& box) const;
    bool inSameMolecule(int atom1, int atom2, const model::MCState& state) const;
    
    // Data members
    std::vector<DrudeParticle> m_particles;
    std::vector<ScreenedPair> m_screenedPairs;
    DrudeSCFParams m_params;
    DrudeAlgorithm m_algorithm;
    
    // Algorithm implementations
    std::unique_ptr<DrudeSCF> m_scfOptimizer;
    std::unique_ptr<DrudeOPT3> m_opt3Optimizer;
    std::unique_ptr<DrudeFBP> m_fbpOptimizer;
    std::unique_ptr<DrudeLBFGS> m_lbfgsOptimizer;
    DrudeOptimizer* m_currentOptimizer;
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc