#ifndef UNIFIED_ACCEPTANCE_HPP
#define UNIFIED_ACCEPTANCE_HPP

#include "CavityBiasCore.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Unified acceptance probability calculator following Metropolis-Hastings
class UnifiedAcceptance {
public:
    // Insertion acceptance probability with proper cavity bias
    // Based on: A_ins = min(1, exp(βμ - βΔE) · V_cav_before / ((n+1)·Λ³))
    static double insertionProbability(
        int n_before,           // Number of molecules before insertion
        double deltaE,          // Energy change (E_after - E_before)
        double beta,            // 1/kT
        double mu,              // Chemical potential
        double V_cavity_before, // Cavity volume BEFORE insertion (nm³)
        double lambda3 = 1.0    // Thermal wavelength cubed (nm³), default 1.0
    ) {
        // Use log-space for numerical stability
        double logProb = beta * mu - beta * deltaE
                       + std::log(V_cavity_before)
                       - std::log(static_cast<double>(n_before + 1))
                       - std::log(lambda3);
        
        return std::min(1.0, std::exp(logProb));
    }
    
    // Deletion acceptance probability with proper cavity bias
    // Based on: A_del = min(1, exp(-βμ - βΔE) · n·Λ³ / V_cav_after)
    static double deletionProbability(
        int n_before,          // Number of molecules before deletion
        double deltaE,         // Energy change (E_after - E_before)
        double beta,           // 1/kT
        double mu,             // Chemical potential
        double V_cavity_after, // Cavity volume AFTER deletion (nm³)
        double lambda3 = 1.0   // Thermal wavelength cubed (nm³)
    ) {
        if (n_before <= 0) return 0.0;
        
        double logProb = -beta * mu - beta * deltaE
                       + std::log(static_cast<double>(n_before))
                       + std::log(lambda3)
                       - std::log(V_cavity_after);
        
        return std::min(1.0, std::exp(logProb));
    }
    
    // Standard GCMC without cavity bias (for comparison)
    static double insertionProbabilityStandard(
        int n_before,
        double deltaE,
        double beta,
        double mu,
        double V_total,        // Total box volume (nm³)
        double lambda3 = 1.0
    ) {
        return insertionProbability(n_before, deltaE, beta, mu, V_total, lambda3);
    }
    
    static double deletionProbabilityStandard(
        int n_before,
        double deltaE,
        double beta,
        double mu,
        double V_total,
        double lambda3 = 1.0
    ) {
        return deletionProbability(n_before, deltaE, beta, mu, V_total, lambda3);
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // UNIFIED_ACCEPTANCE_HPP