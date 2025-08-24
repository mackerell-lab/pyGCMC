#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_MAIN_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_MAIN_HPP

#include <memory>
#include <map>
#include <vector>
#include "../common/MovementResult.hpp"
#include "../common/MovementParams.hpp"
#include "../common/MovementStatistics.hpp"
#include "../common/MovementUtils.hpp"  // For Vector3

namespace pygcmc {

// Forward declarations
namespace model {
namespace montecarlo {
class MCState;
}
}

namespace platform {
namespace cpu {
namespace movement {

// Forward declarations
class ActivePool;
class CavityManager;
class ConfigBiasManager;
class EnergyInterface;

/**
 * Main Movement Module for GCMC operations
 * 
 * This class provides the high-level interface for all GCMC movement operations
 * including insertion, deletion, translation, and rotation moves.
 */
class MovementModule {
public:
    // Constructor and destructor
    MovementModule();
    explicit MovementModule(const MovementParams& params);
    ~MovementModule();
    
    // Disable copy
    MovementModule(const MovementModule&) = delete;
    MovementModule& operator=(const MovementModule&) = delete;
    
    // Main movement functions
    MovementResult attemptInsertion(model::montecarlo::MCState& state, int moleculeType = 0);
    MovementResult attemptDeletion(model::montecarlo::MCState& state, int residueIndex = -1);
    MovementResult attemptTranslation(model::montecarlo::MCState& state, int residueIndex = -1);
    MovementResult attemptRotation(model::montecarlo::MCState& state, int residueIndex = -1);
    
    // Advanced moves with biasing
    MovementResult attemptCavityBiasInsertion(model::montecarlo::MCState& state, int moleculeType = 0);
    MovementResult attemptConfigBiasRotation(model::montecarlo::MCState& state, int residueIndex = -1);
    
    // Multi-insertion CBMC
    std::vector<MovementResult> attemptMultiInsertionCBMC(model::montecarlo::MCState& state, int moleculeType = 0);
    
    // Cavity analysis
    std::vector<Vector3> findCavities(const model::montecarlo::MCState& state);
    double calculateCavityVolume(const model::montecarlo::MCState& state);
    
    // Parameter management
    void setParams(const MovementParams& params);
    MovementParams getParams() const;
    
    // Statistics
    using Statistics = MovementStatistics;
    double calculateAcceptanceRate(const std::string& moveType) const;
    std::map<std::string, Statistics> getStatistics() const;
    void resetStatistics();
    
    // P2: Enhanced statistics access (returns map for binding conversion)
    std::map<std::string, double> getProposalStatsMap() const;
    std::map<std::string, double> getCavityStatsMap() const;
    
    // Active pool management
    ActivePool* getActivePool() { return activePool_.get(); }
    const ActivePool* getActivePool() const { return activePool_.get(); }
    
private:
    // Implementation details
    class Impl;
    std::unique_ptr<Impl> pImpl_;
    
    // Core components
    MovementParams params_;
    std::unique_ptr<ActivePool> activePool_;
    std::unique_ptr<CavityManager> cavityManager_;
    std::unique_ptr<ConfigBiasManager> configBiasManager_;
    std::unique_ptr<EnergyInterface> energyCalc_;
    
#ifdef PYGCMC_USE_PROPOSAL_LAYER
    // Proposal layer (when enabled)
    std::unique_ptr<class ProposalMain> proposalMain_;
#endif
    
    // Statistics tracking
    std::map<std::string, Statistics> stats_;
    
    // Track last inserted residue for paired deletion
    int lastInsertedResidueIndex_ = -1;
    
    // Helper functions
    void initializeComponents();
    void updateStatistics(const std::string& moveType, bool accepted, double energyChange);
    void fillBasicProposalStats(std::map<std::string, double>& result) const;
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_MAIN_HPP