#pragma once

/**
 * @file DrudeStructures.hpp
 * @brief Data structures for Drude oscillator model
 */

#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Physical constants used in Drude calculations
 */
namespace DrudeConstants {
    constexpr double ONE_4PI_EPS0 = 138.935456;  // kJ·nm/mol/e^2
    constexpr double DRUDE_MASS = 0.4;           // amu
}

/**
 * @brief Drude oscillator parameters for a single particle
 * 
 * Following OpenMM's parameter convention for compatibility
 */
struct DrudeParticle {
    int drudeIndex;           // Index of the Drude particle
    int parentIndex;          // Index of the parent atom
    int aniso1Index = -1;     // Index for anisotropy axis 1 (-1 if isotropic)
    int aniso2Index = -1;     // Index for anisotropy axis 2
    int aniso3Index = -1;     // Index for anisotropy axis 3
    int aniso4Index = -1;     // Index for anisotropy axis 4
    double charge;            // Charge on the Drude particle (typically negative)
    double polarizability;    // Isotropic polarizability (nm^3)
    double aniso12 = 1.0;     // Anisotropy scale factor for axis 1-2
    double aniso34 = 1.0;     // Anisotropy scale factor for axis 3-4
    
    // Derived quantities (computed once for efficiency)
    double kSpring;           // Spring constant = charge^2 / (4πε₀α)
    double kAniso1 = 0.0;     // Anisotropic spring constant 1
    double kAniso2 = 0.0;     // Anisotropic spring constant 2
    
    /**
     * @brief Compute derived spring constants
     */
    void computeSpringConstants() {
        // Isotropic spring constant
        // From first principles and OpenMM implementation:
        // The induced dipole moment μ = α × E
        // For Drude model: μ = -q_drude × d (where d is displacement)
        // Force balance: q_drude × E = k × d
        // Therefore: d = q_drude × E / k
        // Combining: α × E = -q_drude × (q_drude × E / k) = q_drude² × E / k
        // Thus: k = q_drude² / α
        // In MD units: k[kJ/mol/nm²] = q²[e²] × ONE_4PI_EPS0[kJ·nm/mol/e²] / α[nm³]
        kSpring = charge * charge * DrudeConstants::ONE_4PI_EPS0 / polarizability;
        
        // Anisotropic contributions (if needed)
        if (aniso1Index >= 0 && aniso2Index >= 0) {
            kAniso1 = kSpring * (aniso12 - 1.0);
        }
        if (aniso3Index >= 0 && aniso4Index >= 0) {
            kAniso2 = kSpring * (aniso34 - 1.0);
        }
    }
};

/**
 * @brief Thole-screened dipole-dipole interaction
 */
struct ScreenedPair {
    int dipole1;      // Index in DrudeParticle array (not atom index)
    int dipole2;      // Index in DrudeParticle array
    double thole;     // Thole screening parameter (typically ~1.3)
};

/**
 * @brief SCF convergence parameters
 * 
 * Default values match OpenMM for consistency
 */
struct DrudeSCFParams {
    double tolerance = 10.0;          // Force tolerance (kJ/mol/nm) - tighter for better convergence
    int maxIterations = 100;          // Maximum SCF iterations
    double dampingFactor = 0.5;       // Damping for stability
    double maxDrudeDistance = 0.02;   // Maximum Drude-parent distance (nm)
};

/**
 * @brief Available optimization algorithms
 */
enum class DrudeAlgorithm {
    SCF,      // Self-Consistent Field iteration
    OPT3,     // 3rd order perturbation theory
    FBP       // Force Balance Predictor
};

/**
 * @brief OPT3 expansion coefficients
 * 
 * For the expansion: r = c0*r0 + c1*r1 + c2*r2 + c3*r3
 * Default is pure first-order response
 */
struct OPT3Coefficients {
    double c0 = 0.0;  // Zero-order (static field only)
    double c1 = 1.0;  // First-order (standard linear response)
    double c2 = 0.0;  // Second-order correction
    double c3 = 0.0;  // Third-order correction
};

/**
 * @brief Compute Thole screening function
 * @param r Distance between dipoles
 * @param alpha_i Polarizability of dipole i
 * @param alpha_j Polarizability of dipole j
 * @param thole Thole parameter
 * @return Screening factor (0 to 1)
 */
inline double computeTholeScreening(double r, double alpha_i, double alpha_j, double thole) {
    // Special case: no screening if thole = 0
    if (thole == 0.0) {
        return 1.0;
    }
    
    // Calculate effective polarizability
    double alpha_eff = std::pow(alpha_i * alpha_j, 1.0/6.0);
    
    // Calculate screening parameter u = thole * r / alpha_eff
    double u = thole * r / alpha_eff;
    
    // Avoid numerical issues for very large u
    if (u > 50.0) {
        return 1.0;
    }
    
    // Calculate screening function S(u) = 1 - (1 + u/2) * exp(-u)
    return 1.0 - (1.0 + u / 2.0) * std::exp(-u);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc