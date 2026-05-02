#pragma once

/**
 * @file DrudeLBFGS.hpp
 * @brief L-BFGS optimizer for Drude oscillators
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include <vector>
#include <deque>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief L-BFGS optimizer for Drude positions
 *
 * Implements Limited-memory BFGS optimization algorithm
 * for finding equilibrium Drude positions. This matches
 * OpenMM's DrudeSCFIntegrator approach for higher precision.
 */
class DrudeLBFGS : public DrudeOptimizer {
public:
    DrudeLBFGS() = default;
    ~DrudeLBFGS() = default;

    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;

    const char* getName() const override { return "L-BFGS"; }

private:
    // L-BFGS memory parameter (typically 5-10)
    static constexpr int MEMORY_SIZE = 5;

    // Line search parameters
    static constexpr double ARMIJO_C1 = 1e-4;
    static constexpr double WOLFE_C2 = 0.9;
    static constexpr int MAX_LINE_SEARCH = 20;

    /**
     * @brief State vector for all Drude positions
     */
    struct StateVector {
        std::vector<double> x;  // Flattened 3N vector

        StateVector() = default;
        StateVector(size_t n) : x(3 * n, 0.0) {}

        void set_size(size_t n) { x.resize(3 * n, 0.0); }
        size_t size() const { return x.size(); }

        double& operator[](size_t i) { return x[i]; }
        const double& operator[](size_t i) const { return x[i]; }

        // Vector operations
        void zero() { std::fill(x.begin(), x.end(), 0.0); }
        double dot(const StateVector& other) const;
        void axpy(double a, const StateVector& y);  // x = x + a*y
        void scale(double a);  // x = a*x
        double norm() const;
    };

    /**
     * @brief History entry for L-BFGS
     */
    struct HistoryEntry {
        StateVector s;  // x_{k+1} - x_k
        StateVector y;  // g_{k+1} - g_k
        double rho;     // 1 / (y^T s)
    };

    // Helper functions
    void packState(const model::MCState& state,
                   const std::vector<DrudeParticle>& particles,
                   StateVector& x) const;

    void unpackState(const StateVector& x,
                     const std::vector<DrudeParticle>& particles,
                     model::MCState& state) const;

    double evaluateEnergy(const StateVector& x,
                         const std::vector<DrudeParticle>& particles,
                         const std::vector<ScreenedPair>& screenedPairs,
                         model::MCState& state) const;

    void evaluateGradient(const StateVector& x,
                         const std::vector<DrudeParticle>& particles,
                         const std::vector<ScreenedPair>& screenedPairs,
                         model::MCState& state,
                         StateVector& gradient) const;

    void computeSearchDirection(const StateVector& gradient,
                               const std::deque<HistoryEntry>& history,
                               double H0,
                               StateVector& direction) const;

    double lineSearch(StateVector& x,
                     const StateVector& direction,
                     double f0,
                     const StateVector& g0,
                     const std::vector<DrudeParticle>& particles,
                     const std::vector<ScreenedPair>& screenedPairs,
                     model::MCState& state,
                     StateVector& gradient) const;

    // Force calculation (reuse from SCF)
    void calculateElectricField(const model::MCState& state,
                               const std::vector<DrudeParticle>& particles,
                               const std::vector<ScreenedPair>& screenedPairs,
                               std::vector<Vec3>& electricField) const;

    void calculateExternalField(const model::MCState& state,
                               const std::vector<DrudeParticle>& particles,
                               std::vector<Vec3>& electricField) const;

    void calculateInducedField(const model::MCState& state,
                              const std::vector<DrudeParticle>& particles,
                              const std::vector<ScreenedPair>& screenedPairs,
                              std::vector<Vec3>& electricField) const;

    void applyPBC(double& dx, double& dy, double& dz,
                  const std::array<double, 3>& box) const;

    bool inSameMolecule(int atom1, int atom2,
                       const model::MCState& state) const;
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
