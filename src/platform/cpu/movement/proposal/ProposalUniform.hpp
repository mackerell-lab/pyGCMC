// ============================================================================
// ProposalUniform.hpp - Uniform proposal sampler
// ============================================================================

#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_UNIFORM_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_UNIFORM_HPP

#include "ProposalInterface.hpp"
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Uniform proposal sampler
 * Samples uniformly within the simulation box
 */
class ProposalUniform : public ProposalInterface {
private:
    mutable std::mt19937 rng_;
    mutable std::uniform_real_distribution<double> dist_;
    
public:
    ProposalUniform(unsigned seed = std::random_device{}());
    
    Vector3 propose(const MCState& state, ProposalInfo& info) override;
    
    ProposalType getType() const override { 
        return ProposalType::Uniform; 
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_UNIFORM_HPP