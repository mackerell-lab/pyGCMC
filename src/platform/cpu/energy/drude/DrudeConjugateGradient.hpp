#ifndef DRUDE_CONJUGATE_GRADIENT_HPP
#define DRUDE_CONJUGATE_GRADIENT_HPP

#include "DrudeForce.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

// Forward declaration
class DrudeForce;

/**
 * Conjugate Gradient solver for Drude oscillator positions
 * 
 * Based on the linear response theory, the equilibrium Drude positions
 * satisfy the linear system:
 *   K * delta_r = F_external
 * 
 * where:
 * - K is the effective force constant matrix (includes spring + Thole)
 * - delta_r is the displacement of Drude from parent
 * - F_external is the external electric field force
 * 
 * This implementation uses a matrix-free conjugate gradient method
 * that avoids explicitly constructing the full matrix K.
 */
class DrudeConjugateGradient {
public:
    /**
     * Minimize Drude positions using Conjugate Gradient method
     * 
     * @param state MCState containing atomic positions
     * @param drudeForce DrudeForce object with particle definitions
     * @param tolerance Convergence tolerance in kJ/mol/nm
     * @param maxIterations Maximum CG iterations
     * @return true if converged, false otherwise
     */
    static bool minimizeDrudePositions(
        model::MCState& state,
        const DrudeForce& drudeForce,
        double tolerance = 1.0,
        int maxIterations = 100
    );

private:
    /**
     * Apply the system matrix K to a vector x
     * K*x includes spring forces and Thole interactions
     * 
     * This is the key operation for matrix-free CG
     */
    static std::vector<DrudeForce::Vec3> applySystemMatrix(
        const model::MCState& state,
        const DrudeForce& drudeForce,
        const std::vector<DrudeForce::Vec3>& x
    );
    
    /**
     * Calculate the right-hand side vector b = F_external
     * This is the negative of the total force on Drude particles
     * (excluding spring and Thole forces between Drudes)
     */
    static std::vector<DrudeForce::Vec3> calculateRHS(
        const model::MCState& state,
        const DrudeForce& drudeForce
    );
    
    /**
     * Compute dot product of two vector arrays
     */
    static double dotProduct(
        const std::vector<DrudeForce::Vec3>& a,
        const std::vector<DrudeForce::Vec3>& b
    );
    
    /**
     * Update Drude positions based on solution vector
     */
    static void updatePositions(
        model::MCState& state,
        const DrudeForce& drudeForce,
        const std::vector<DrudeForce::Vec3>& solution
    );
    
    /**
     * Calculate maximum force magnitude for convergence check
     */
    static double calculateMaxForce(
        const model::MCState& state,
        const DrudeForce& drudeForce
    );
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // DRUDE_CONJUGATE_GRADIENT_HPP