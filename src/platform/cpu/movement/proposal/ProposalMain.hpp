// ============================================================================
// ProposalMain.hpp - Main manager for proposal system
// ============================================================================

#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_MAIN_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_MAIN_HPP

#include "ProposalInterface.hpp"
#include "ProposalTypes.hpp"
#include "../common/MovementParams.hpp"
#include "../bias/CavityBias.hpp"
#include <memory>

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
    
    // Statistics
    struct Statistics {
        int totalAttempts = 0;
        int acceptedMoves = 0;
        double acceptanceRate = 0.0;
        
        void update(bool accepted) {
            totalAttempts++;
            if (accepted) acceptedMoves++;
            acceptanceRate = (totalAttempts > 0) ? 
                           static_cast<double>(acceptedMoves) / totalAttempts : 0.0;
        }
        
        void reset() {
            totalAttempts = 0;
            acceptedMoves = 0;
            acceptanceRate = 0.0;
        }
    };
    Statistics stats_;
    
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
     * Get statistics
     */
    const Statistics& getStatistics() const { return stats_; }
    
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