// ============================================================================
// ProposalInterface.hpp - Base interface for all proposal strategies
// ============================================================================

#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_INTERFACE_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_INTERFACE_HPP

#include "ProposalTypes.hpp"
#include "../common/MovementUtils.hpp"

// Forward declarations
namespace pygcmc {
namespace model {
namespace montecarlo {
    struct MCState;
}
}
}

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCState;

/**
 * Base interface for all proposal strategies
 */
class ProposalInterface {
public:
    virtual ~ProposalInterface() = default;

    /**
     * Generate a proposal position
     * @param state Current MC state
     * @param info Output proposal information
     * @return Proposed position in nm
     */
    virtual Vector3 propose(const MCState& state, ProposalInfo& info) = 0;

    /**
     * Update internal state after move
     * @param accepted Whether the move was accepted
     * @param info The proposal that was tested
     */
    virtual void update(bool accepted, const ProposalInfo& info) {}

    /**
     * Reset internal statistics
     */
    virtual void reset() {}

    /**
     * Get proposal type
     */
    virtual ProposalType getType() const = 0;
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_INTERFACE_HPP
