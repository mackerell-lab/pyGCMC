#ifndef PYGCMC_IO_OUTPUT_TRAJECTORY_WRITER_HPP
#define PYGCMC_IO_OUTPUT_TRAJECTORY_WRITER_HPP

#include "../../model/montecarlo/MCMain.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace pygcmc {
namespace io {
namespace output {

/**
 * @brief Trajectory writer for simulation output
 * 
 * This class handles writing trajectory files in various formats including
 * PDB, XYZ, and DAT formats. It supports both single-frame and multi-frame
 * outputs with proper element inference and formatting.
 */
class TrajectoryWriter {
public:
    /**
     * @brief Configuration for trajectory writing
     */
    struct Config {
        std::string format = "pdb";        // Output format: "pdb", "xyz", "dat"
        std::string prefix = "gcmc";       // Filename prefix
        bool compressOutput = false;       // Whether to compress output (future)
        int precision = 3;                 // Decimal precision for coordinates
        bool multiFrame = false;           // Write multiple frames to same file
        bool continuousFile = false;       // Keep file open between writes
        
        Config() = default;
    };
    
    // Constructor and destructor
    explicit TrajectoryWriter(const Config& config);
    ~TrajectoryWriter();
    
    // File operations
    bool open(const std::string& filename = "");
    void close();
    bool isOpen() const { return isOpen_; }
    
    // Main write function
    bool writeTrajectory(const model::montecarlo::MCState& state,
                        int step,
                        const std::string& filename = "");
    
    // Configuration
    void setConfig(const Config& config) { config_ = config; }
    const Config& getConfig() const { return config_; }
    
    // Frame counting
    int getFrameCount() const { return frameCount_; }
    void resetFrameCount() { frameCount_ = 0; }
    
private:
    Config config_;
    std::ofstream file_;
    std::string currentFilename_;
    int frameCount_;
    bool isOpen_;
    
    // Format-specific writers
    bool writePDB(const model::montecarlo::MCState& state,
                 int step,
                 std::ofstream& out);
    
    bool writeXYZ(const model::montecarlo::MCState& state,
                 int step,
                 std::ofstream& out);
    
    bool writeDAT(const model::montecarlo::MCState& state,
                 int step,
                 std::ofstream& out);
    
    bool writeTOP(const model::montecarlo::MCState& state,
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
                             double x, double y, double z,
                             double occupancy = 1.0,
                             double tempFactor = 0.0,
                             const std::string& element = "") const;
};

/**
 * @brief Data writer for analysis output
 * 
 * This class provides simple columnar data output for analysis,
 * supporting headers and numerical data rows.
 */
class DataWriter {
public:
    explicit DataWriter(const std::string& filename);
    ~DataWriter();
    
    void close();
    void writeHeader(const std::vector<std::string>& columns);
    void writeRow(const std::vector<double>& values);
    void writeComment(const std::string& comment);
    
private:
    std::ofstream file_;
    bool headerWritten_;
};

} // namespace output
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_OUTPUT_TRAJECTORY_WRITER_HPP