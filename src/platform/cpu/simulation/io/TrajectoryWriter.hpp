#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_IO_TRAJECTORY_WRITER_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_IO_TRAJECTORY_WRITER_HPP

#include "../../../../model/montecarlo/MCMain.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

/**
 * @brief Trajectory writer for simulation output
 * 
 * This class handles writing trajectory files in various formats.
 */
class TrajectoryWriter {
public:
    /**
     * @brief Configuration for trajectory writing
     */
    struct Config {
        std::string format;
        std::string prefix;
        bool compressOutput;
        int precision;
        
        Config() : format("pdb"), prefix("gcmc"), compressOutput(false), precision(3) {}
    };
    
    // Constructor
    TrajectoryWriter(const Config& config);
    ~TrajectoryWriter();
    
    // Write trajectory
    bool writeTrajectory(const model::montecarlo::MCState& state,
                        int step,
                        const std::string& filename = "");
    
    // Configuration
    void setConfig(const Config& config) { config_ = config; }
    const Config& getConfig() const { return config_; }
    
private:
    Config config_;
    
    // Format-specific writers
    bool writePDB(const model::montecarlo::MCState& state,
                 int step,
                 std::ofstream& out);
    
    bool writeXYZ(const model::montecarlo::MCState& state,
                 int step,
                 std::ofstream& out);
    
    // Helper methods
    std::string generateFilename(int step) const;
    void writePDBHeader(std::ofstream& out, int step) const;
    void writePDBBox(std::ofstream& out, 
                    const std::vector<double>& box) const;
    std::string formatPDBAtom(int serial,
                             const std::string& atomName,
                             const std::string& resName,
                             int resSeq,
                             double x, double y, double z) const;
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_IO_TRAJECTORY_WRITER_HPP