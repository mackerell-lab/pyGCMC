// ============================================================================
// ProposalMain.hpp - Main manager for proposal system
// ============================================================================

#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_MAIN_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_MAIN_HPP

#include "ProposalInterface.hpp"
#include "ProposalTypes.hpp"
#include "ProposalStatistics.hpp"
#include "../common/MovementParams.hpp"
#include "../bias/CavityBias.hpp"
#include <memory>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Main proposal system manager
 */
class ProposalMain {
private:
    std::unique_ptr<ProposalInterface> sampler_;
    ProposalType currentType_;
    const MovementParams& params_;
    std::shared_ptr<CavityManager> cavityManager_;

    // P2: Enhanced statistics
    ProposalStatistics stats_;

    // Timing helpers
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;

public:
    ProposalMain(const MovementParams& params,
                 std::shared_ptr<CavityManager> cavityManager = nullptr);

    /**
     * Generate a proposal
     */
    ProposalInfo generateProposal(const MCState& state);

    /**
     * Update after move
     */
    void updateState(bool accepted, const ProposalInfo& info);

    /**
     * Switch proposal type
     */
    void switchType(ProposalType type);

    /**
     * Get current type
     */
    ProposalType getCurrentType() const { return currentType_; }

    /**
     * Get statistics (P2 enhanced)
     */
    const ProposalStatistics& getStatistics() const { return stats_; }

    /**
     * Get cavity manager statistics (if available)
     */
    CavityManager::Statistics getCavityStatistics() const {
        if (cavityManager_) {
            return cavityManager_->getStatistics();
        }
        return CavityManager::Statistics();
    }

    /**
     * Check if should switch mode (adaptive)
     */
    bool shouldSwitchMode(const MCState& state) const;

    /**
     * Reset all
     */
    void reset();

private:
    ProposalType determineType(const MovementParams& params);
    std::unique_ptr<ProposalInterface> createSampler(ProposalType type);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_MAIN_HPP
