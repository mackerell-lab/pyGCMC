#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_IO_CHECKPOINT_MANAGER_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_IO_CHECKPOINT_MANAGER_HPP

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../stats/StatisticsTracker.hpp"
#include <string>
#include <fstream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

/**
 * @brief Checkpoint manager for simulation save/load
 * 
 * This class handles saving and loading simulation checkpoints.
 */
class CheckpointManager {
public:
    /**
     * @brief Configuration for checkpoint operations
     */
    struct Config {
        std::string prefix;
        bool compress;
        int version;
        
        Config() : prefix("gcmc"), compress(false), version(1) {}
    };
    
    // Constructor
    CheckpointManager(const Config& config);
    ~CheckpointManager();
    
    // Save checkpoint
    bool saveCheckpoint(const model::montecarlo::MCState& state,
                       const StatisticsTracker& stats,
                       int step,
                       const std::string& filename = "");
    
    // Load checkpoint
    bool loadCheckpoint(model::montecarlo::MCState& state,
                       StatisticsTracker& stats,
                       int& step,
                       const std::string& filename);
    
    // Validate checkpoint file
    bool validateCheckpoint(const std::string& filename) const;
    
    // Configuration
    void setConfig(const Config& config) { config_ = config; }
    const Config& getConfig() const { return config_; }
    
private:
    Config config_;
    
    // Helper methods
    std::string generateFilename(int step) const;
    bool writeHeader(std::ofstream& out, int version) const;
    bool readHeader(std::ifstream& in, int& version) const;
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_IO_CHECKPOINT_MANAGER_HPP