// ============================================================================
// ProposalCavity.hpp - Cavity-biased proposal sampler
// ============================================================================

#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_CAVITY_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_CAVITY_HPP

#include "ProposalInterface.hpp"
#include "ProposalUniform.hpp"
#include "../bias/CavityBias.hpp"
#include <memory>
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Cavity-biased proposal sampler
 */
class ProposalCavity : public ProposalInterface {
private:
    std::shared_ptr<CavityManager> cavityManager_;
    std::unique_ptr<ProposalUniform> uniformFallback_;
    mutable std::mt19937 rng_;
    mutable std::uniform_real_distribution<double> dist_;
    
public:
    ProposalCavity(std::shared_ptr<CavityManager> cavityManager, 
                   unsigned seed = std::random_device{}());
    
    Vector3 propose(const MCState& state, ProposalInfo& info) override;
    
    ProposalType getType() const override { 
        return ProposalType::Cavity; 
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_CAVITY_HPP