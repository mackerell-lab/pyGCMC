#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_COMMON_GRANDCANONICALTERMS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_COMMON_GRANDCANONICALTERMS_HPP

#include <limits>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Log-form representation of all factors entering the grand canonical
 *        acceptance ratio for insertion/deletion moves.
 */
struct GrandCanonicalTerms {
    int speciesId = -1;
    int countBefore = 0;   // Particle count prior to the move
    int countAfter = 0;    // Particle count after the move

    double beta = 0.0;               // 1/(k_B T)
    double chemicalPotential = 0.0;  // μ in kJ/mol
    double deltaEnergy = 0.0;        // E_after - E_before in kJ/mol

    double logVolume = 0.0;          // log(V) with V in nm^3
    double logLambda3 = 0.0;         // log(Λ^3) with Λ in nm

    double logProposalForward = 0.0; // log(q_forward)
    double logProposalReverse = 0.0; // log(q_reverse)

    double logCavityForward = 0.0;   // log(cavity factor forward)
    double logCavityReverse = 0.0;   // log(cavity factor reverse)

    double logRosenbluthForward = 0.0; // log(W_forward/K_forward)
    double logRosenbluthReverse = 0.0; // log(W_reverse/K_reverse)

    double logExtraForward = 0.0;    // Additional bias terms (forward)
    double logExtraReverse = 0.0;    // Additional bias terms (reverse)

    void reset() { *this = GrandCanonicalTerms(); }
};

/**
 * @brief Result of evaluating the acceptance ratio in log-space.
 */
struct GrandCanonicalEvaluation {
    double probability = 0.0;
    double logRatio = -std::numeric_limits<double>::infinity();
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_COMMON_GRANDCANONICALTERMS_HPP
