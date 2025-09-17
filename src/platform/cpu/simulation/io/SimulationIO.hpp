#pragma once

/**
 * @file SimulationIO.hpp
 * @brief I/O operations for GCMC simulation
 */

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../../../io/output/TrajectoryWriter.hpp"
#include "../../../../io/output/CheckpointManager.hpp"
#include "../../../../io/parameters/InpParserMain.hpp"
#include "../stats/StatisticsTracker.hpp"
#include <memory>
#include <string>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace io {

/**
 * @brief Handles all I/O operations for the simulation
 */
class SimulationIO {
public:
    /**
     * @brief Configuration for I/O operations
     */
    struct Config {
        std::string outputPrefix = "gcmc";
        std::string trajectoryFormat = "pdb";  // pdb, xyz, dat
        bool compressTrajectory = false;
        bool writeEnergies = true;
        bool writeDensities = true;
        int precision = 6;
        bool appendMode = false;
        
        Config() = default;
    };
    
    // Constructor and destructor
    explicit SimulationIO(const Config& config);
    ~SimulationIO();
    
    // Input operations
    bool readInputFile(const std::string& filename,
                      model::montecarlo::MCState& state,
                      model::param::Param& params);
    
    bool readFragmentFile(const std::string& filename,
                         std::vector<model::montecarlo::MCResidue>& fragments);
    
    bool readCheckpoint(const std::string& filename,
                       model::montecarlo::MCState& state,
                       int& step);
    
    // Output operations
    void writeTrajectory(const model::montecarlo::MCState& state,
                        int step);
    
    void writeCheckpoint(const model::montecarlo::MCState& state,
                        int step,
                        const StatisticsTracker* stats = nullptr);
    
    void writeStatistics(const StatisticsTracker& stats,
                        int step);
    
    void writeEnergies(const model::montecarlo::MCState& state,
                      int step);
    
    void writeDensities(const std::map<std::string, double>& densities,
                       int step);
    
    void writeFinalReport(const model::montecarlo::MCState& state,
                         const StatisticsTracker& stats);
    
    // File management
    void openTrajectoryFile(const std::string& filename = "");
    void closeTrajectoryFile();
    bool isTrajectoryOpen() const;
    
    void openStatisticsFile(const std::string& filename = "");
    void closeStatisticsFile();
    
    void openEnergyFile(const std::string& filename = "");
    void closeEnergyFile();
    
    // Configuration
    const Config& getConfig() const { return config_; }
    void updateConfig(const Config& config);
    
    // Utilities
    std::string generateFilename(const std::string& suffix,
                                const std::string& extension) const;
    
    static std::string formatTimestamp();
    
private:
    Config config_;
    
    std::unique_ptr<pygcmc::io::output::TrajectoryWriter> trajectoryWriter_;
    std::unique_ptr<pygcmc::io::output::CheckpointManager> checkpointManager_;
    std::unique_ptr<pygcmc::io::output::DataWriter> statisticsWriter_;
    std::unique_ptr<pygcmc::io::output::DataWriter> energyWriter_;
    std::unique_ptr<pygcmc::io::output::DataWriter> densityWriter_;
    
    bool trajectoryOpen_;
    bool statisticsOpen_;
    bool energyOpen_;
    
    // Helper methods
    void ensureDirectoryExists(const std::string& path) const;
    void writeHeader(pygcmc::io::output::DataWriter* writer,
                    const std::vector<std::string>& columns);
};

} // namespace io
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc