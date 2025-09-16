#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_MODULAR_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_MODULAR_HPP

#include "core/SimulationRunner.hpp"
#include "setup/SystemInitializer.hpp"
#include "io/TrajectoryWriter.hpp"
#include "io/CheckpointManager.hpp"
#include "stats/StatisticsTracker.hpp"
#include "../movement/gcmc/GCMCEngine.hpp"
#include "../movement/gcmc/GCMCAcceptance.hpp"
#include "../movement/reservoir/MultiTypeReservoir.hpp"
#include "../../../model/param/ParamMain.hpp"
#include "../../../model/montecarlo/MCMain.hpp"
#include <memory>
#include <string>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

/**
 * @brief Modular GCMC simulation class
 * 
 * This is a refactored version that uses separate modules for:
 * - Core simulation logic (SimulationRunner)
 * - System setup (SystemInitializer)
 * - I/O operations (TrajectoryWriter, CheckpointManager)
 * - Statistics tracking (StatisticsTracker)
 */
class GCMCSimulationModular {
public:
    /**
     * @brief Configuration for the simulation
     */
    struct Config {
        std::string inputFile;
        std::string outputPrefix = "gcmc";
        int printFrequency = 1000;
        int trajectoryFrequency = 10000;
        int checkpointFrequency = 100000;
        bool verbose = false;
        int randomSeed = -1;
        bool enableStatistics = true;
        int statisticsInterval = 1000;
        double convergenceTolerance = 0.01;
    };
    
    // Constructor and destructor
    explicit GCMCSimulationModular(const Config& config);
    ~GCMCSimulationModular();
    
    // Main simulation methods
    bool initialize();
    bool run();
    void finalize();
    
    // Control methods
    void stop();
    bool isRunning() const;
    
    // Analysis methods
    const StatisticsTracker& getStatistics() const;
    void printStatistics() const;
    
    // Trajectory and checkpoint methods (delegated to IO module)
    void writeTrajectory(int step);
    void writeCheckpoint(int step);
    bool loadCheckpoint(const std::string& filename);
    
    // Configuration access
    const Config& getConfig() const { return config_; }
    void updateConfig(const Config& config);
    
    // Fragment information access
    std::vector<SystemInitializer::FragmentConfig> getFragmentInfo() const { 
        return fragmentConfigs_; 
    }
    
private:
    // Configuration
    Config config_;
    bool initialized_;
    bool running_;
    
    // Modules
    std::unique_ptr<SimulationRunner> core_;
    std::unique_ptr<SystemInitializer> setup_;
    std::unique_ptr<TrajectoryWriter> trajectoryWriter_;
    std::unique_ptr<CheckpointManager> checkpointManager_;
    
    // Core components
    std::unique_ptr<model::param::Param> params_;
    std::unique_ptr<model::montecarlo::MCState> state_;
    std::unique_ptr<movement::gcmc::GCMCEngine> engine_;
    std::unique_ptr<movement::gcmc::GCMCAcceptance> acceptance_;
    std::unique_ptr<movement::MultiTypeReservoir> reservoir_;
    
    // Fragment configuration
    std::vector<SystemInitializer::FragmentConfig> fragmentConfigs_;
    
    // Helper methods
    void log(const std::string& message) const;
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_MODULAR_HPP