// ============================================================================
// ProposalTypes.hpp - Common types for proposal system
// ============================================================================

#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_TYPES_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_TYPES_HPP

#include "../common/MovementUtils.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Proposal strategy type
 */
enum class ProposalType {
    Uniform = 0,    // Uniform sampling
    Cavity = 1,     // Cavity-biased
    Color = 2,      // Color-class weighted
    Cluster = 3,    // Cluster-based
    Adaptive = 4    // Adaptive selection
};

/**
 * Proposal information structure
 */
struct ProposalInfo {
    // Core fields
    Vector3 position;                     // Position in nm
    double normalizationFactor = 1.0;     // q_norm for detailed balance
    ProposalType type = ProposalType::Uniform;

    // Cavity-specific
    bool usedCavity = false;
    double cavityBiasFactor = 1.0;        // Ntotal/Ncav

    // Multi-insertion
    int numProposals = 1;                 // M_proposal
    double regionVolume = -1.0;            // V_region in nm³

    /**
     * Calculate normalization for acceptance
     */
    double calculateNormalization() const {
        if (numProposals > 1 && regionVolume > 0) {
            return static_cast<double>(numProposals) * regionVolume;
        } else {
            return 1.0 / cavityBiasFactor;
        }
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_TYPES_HPP
