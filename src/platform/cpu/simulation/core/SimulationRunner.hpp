#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_CORE_SIMULATION_RUNNER_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_CORE_SIMULATION_RUNNER_HPP

#include "../../movement/gcmc/GCMCEngine.hpp"
#include "../../movement/gcmc/GCMCAcceptance.hpp"
#include "../../movement/reservoir/MultiTypeReservoir.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../stats/StatisticsTracker.hpp"
#include <random>
#include <memory>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

// Forward declarations
class GCMCSimulation;

/**
 * @brief Core simulation execution module
 * 
 * This class handles the main simulation loop and MC step execution.
 * It's responsible for performing moves and managing the simulation state.
 */
class SimulationRunner {
public:
    // Move types
    enum MoveType {
        INSERT,
        DELETE,
        TRANSLATE,
        ROTATE
    };
    
    // Constructor
    SimulationRunner(GCMCSimulation* parent);
    ~SimulationRunner();
    
    // Initialize core components
    void initialize(model::montecarlo::MCState* state,
                   movement::MultiTypeReservoir* reservoir,
                   movement::gcmc::GCMCEngine* engine,
                   movement::gcmc::GCMCAcceptance* acceptance);
    
    // Main simulation methods
    bool runSimulation(int numSteps);
    bool performMCStep();
    
    // Move selection
    MoveType selectMoveType();
    int selectFragmentType();
    int selectActiveFragment();
    
    // Statistics access
    StatisticsTracker& getStatistics() { return statistics_; }
    const StatisticsTracker& getStatistics() const { return statistics_; }
    
    // Configuration
    void setSeed(unsigned int seed);
    void setVerbose(bool verbose) { verbose_ = verbose; }
    
    // Energy calculation
    double calculateSystemEnergy();
    
private:
    // Parent simulation object
    GCMCSimulation* parent_;
    
    // Core components (not owned, just pointers)
    model::montecarlo::MCState* state_;
    movement::MultiTypeReservoir* reservoir_;
    movement::gcmc::GCMCEngine* engine_;
    movement::gcmc::GCMCAcceptance* acceptance_;
    
    // Statistics tracking
    StatisticsTracker statistics_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;
    
    // Configuration
    bool verbose_;
    int printFrequency_;
    int trajectoryFrequency_;
    int checkpointFrequency_;
    bool enableStatistics_;
    int statisticsInterval_;
    
    // Helper methods
    bool performInsertion(int fragType);
    bool performDeletion();
    bool performTranslation();
    bool performRotation();
    void updateEnergyStatistics();
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_CORE_SIMULATION_RUNNER_HPP