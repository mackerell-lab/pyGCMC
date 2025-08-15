/**
 * @file DrudeSCFOpenMMExact.hpp
 * @brief OpenMM-exact Drude SCF implementation
 * 
 * This implementation exactly matches OpenMM's DrudeForce behavior:
 * - Complete force expression with dS/dr terms
 * - Energy minimization instead of force balance
 * - Identical numerical parameters and convergence criteria
 */

#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include "../DrudeStructures.hpp"
#include "model/montecarlo/MCMain.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

/**
 * @brief OpenMM-exact parameters
 * All parameters locked to match OpenMM defaults
 */
struct DrudeSCFOpenMMExactParams {
    // Convergence (matches OpenMM DrudeSCFIntegrator)
    double minimizationErrorTolerance = 1e-8;  // Force tolerance
    double maxDrudeDistance = 0.02;            // 0.02 nm hard limit
    int maxIterations = 100;                   // More iterations for energy min
    
    // Energy minimization parameters
    double stepSize = 0.001;                   // Initial step size
    double energyTolerance = 1e-10;           // Energy convergence
    
    // Numerical precision
    bool useDoublePrecision = true;           // Always double
    bool fixedSummationOrder = true;          // Deterministic
    
    // Physical constants (exact OpenMM values)
    static constexpr double ONE_4PI_EPS0 = 138.935456;  // kJ·nm/(mol·e²)
    
    // No softening or approximations
    bool allowSoftening = false;              // Always false
    bool useS1Only = false;                   // Always include dS/dr
    
    int logLevel = 0;                         // 0=silent, 1=summary, 2=detailed
};

/**
 * @brief S1 screening function and its derivative
 */
struct S1Derivatives {
    double S1;       // S1(u) value
    double dS1_du;   // dS1/du derivative
    
    static S1Derivatives calculate(double u) {
        S1Derivatives result;
        
        if (u > 50.0) {
            result.S1 = 1.0;
            result.dS1_du = 0.0;
            return result;
        }
        
        const double exp_u = std::exp(-u);
        
        // S1(u) = 1 - (1 + u/2) * exp(-u)
        result.S1 = 1.0 - (1.0 + 0.5*u) * exp_u;
        
        // dS1/du = d/du[1 - (1 + u/2)*exp(-u)]
        //        = -0.5*exp(-u) + (1 + u/2)*exp(-u) 
        //        = exp(-u) * [(1 + u/2) - 0.5]
        //        = exp(-u) * (0.5 + u/2)
        //        = (1 + u)/2 * exp(-u)
        result.dS1_du = 0.5 * (1.0 + u) * exp_u;
        
        return result;
    }
};

/**
 * @brief OpenMM-exact Drude SCF optimizer
 * 
 * Implements energy minimization with complete force expressions
 * to exactly match OpenMM's DrudeForce behavior
 */
class DrudeSCFOpenMMExact {
public:
    using Vec3 = std::array<double, 3>;
    
    DrudeSCFOpenMMExact() = default;
    ~DrudeSCFOpenMMExact() = default;
    
    /**
     * @brief Optimize Drude positions using energy minimization
     * @return true if converged within tolerance
     */
    bool optimize(model::montecarlo::MCState& state,
                  const std::vector<DrudeParticle>& particles,
                  const std::vector<ScreenedPair>& pairs,
                  const DrudeSCFOpenMMExactParams& params = {});
    
    /**
     * @brief Calculate total energy (for verification)
     */
    double calculateEnergy(const model::montecarlo::MCState& state,
                          const std::vector<DrudeParticle>& particles,
                          const std::vector<ScreenedPair>& pairs) const;
    
    /**
     * @brief Get iteration count from last optimization
     */
    int getLastIterationCount() const { return m_lastIterationCount; }
    
private:
    /**
     * @brief Calculate total system energy
     * U_total = U_spring + U_coulomb + U_screening
     */
    double calculateTotalEnergy(const model::montecarlo::MCState& state,
                                const std::vector<DrudeParticle>& particles,
                                const std::vector<ScreenedPair>& pairs) const;
    
    /**
     * @brief Calculate forces with complete expression
     * F = -∇U including dS/dr terms
     */
    void calculateForces(const model::montecarlo::MCState& state,
                        const std::vector<DrudeParticle>& particles,
                        const std::vector<ScreenedPair>& pairs,
                        std::vector<Vec3>& forces) const;
    
    /**
     * @brief Calculate screened Coulomb energy between two charges
     */
    double calculateScreenedEnergy(const Vec3& r_i, const Vec3& r_j,
                                   double q_i, double q_j,
                                   double alpha_eff, double thole,
                                   const std::array<double, 3>& box) const;
    
    /**
     * @brief Calculate screened Coulomb force with dS/dr term
     * F = -q_i*q_j*[S1/r² + (dS1/du)(du/dr)/r] * r_hat
     */
    Vec3 calculateScreenedForce(const Vec3& r_i, const Vec3& r_j,
                                double q_i, double q_j,
                                double alpha_eff, double thole,
                                const std::array<double, 3>& box) const;
    
    /**
     * @brief Calculate spring energy for all Drude particles
     */
    double calculateSpringEnergy(const model::montecarlo::MCState& state,
                                 const std::vector<DrudeParticle>& particles) const;
    
    /**
     * @brief Calculate spring forces
     */
    void calculateSpringForces(const model::montecarlo::MCState& state,
                               const std::vector<DrudeParticle>& particles,
                               std::vector<Vec3>& forces) const;
    
    /**
     * @brief Calculate unscreened Coulomb contributions
     */
    double calculateUnscreenedEnergy(const model::montecarlo::MCState& state,
                                     const std::vector<DrudeParticle>& particles) const;
    
    /**
     * @brief Apply periodic boundary conditions
     */
    void applyPBC(double& dx, double& dy, double& dz,
                  const std::array<double, 3>& box) const;
    
    /**
     * @brief Check if two atoms are in the same molecule
     */
    bool inSameMolecule(int atom1, int atom2, const model::montecarlo::MCState& state) const;
    
    /**
     * @brief Update Drude positions using gradient descent
     */
    bool updatePositions(model::montecarlo::MCState& state,
                        const std::vector<DrudeParticle>& particles,
                        const std::vector<Vec3>& forces,
                        double stepSize,
                        const DrudeSCFOpenMMExactParams& params) const;
    
    /**
     * @brief Line search for optimal step size
     */
    double lineSearch(const model::montecarlo::MCState& state,
                     const std::vector<DrudeParticle>& particles,
                     const std::vector<ScreenedPair>& pairs,
                     const std::vector<Vec3>& searchDirection,
                     double initialStep) const;
    
    /**
     * @brief L-BFGS optimization (advanced)
     */
    bool optimizeLBFGS(model::montecarlo::MCState& state,
                      const std::vector<DrudeParticle>& particles,
                      const std::vector<ScreenedPair>& pairs,
                      const DrudeSCFOpenMMExactParams& params);
    
    // State tracking
    mutable int m_lastIterationCount = 0;
    mutable double m_lastEnergy = 0.0;
    mutable double m_lastForceNorm = 0.0;
    
    // L-BFGS history (if used)
    struct LBFGSHistory {
        std::vector<std::vector<double>> s;  // Position differences
        std::vector<std::vector<double>> y;  // Gradient differences
        std::vector<double> rho;             // 1/(y·s)
        int m = 10;                          // History size
    };
    std::unique_ptr<LBFGSHistory> m_lbfgsHistory;
};

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc