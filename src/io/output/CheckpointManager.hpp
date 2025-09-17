#ifndef PYGCMC_IO_OUTPUT_CHECKPOINT_MANAGER_HPP
#define PYGCMC_IO_OUTPUT_CHECKPOINT_MANAGER_HPP

#include "../../model/montecarlo/MCMain.hpp"
#include "../../platform/cpu/simulation/stats/StatisticsTracker.hpp"
#include <string>
#include <fstream>

namespace pygcmc {
namespace io {
namespace output {

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
                       const platform::cpu::simulation::StatisticsTracker& stats,
                       int step,
                       const std::string& filename = "");
    
    // Load checkpoint
    bool loadCheckpoint(model::montecarlo::MCState& state,
                       platform::cpu::simulation::StatisticsTracker& stats,
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

} // namespace output
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_OUTPUT_CHECKPOINT_MANAGER_HPP