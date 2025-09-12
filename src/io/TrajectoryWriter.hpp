// src/io/TrajectoryWriter.hpp
#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace io {

class TrajectoryWriter {
public:
    enum class Format {
        PDB,
        XYZ,
        DAT,  // Simple data format for analysis
        TOP   // Topology format
    };

    TrajectoryWriter(const std::string& filename, Format format = Format::PDB);
    ~TrajectoryWriter();

    // Write a single frame
    void writeFrame(const model::MCState& state, int frameNumber = -1);
    
    // Write statistics data
    void writeStatistics(const model::MCInfo::Statistics& stats);
    
    // Write topology (for TOP format)
    void writeTopology(const model::MCState& state);
    
    // Close file
    void close();

private:
    void writePDBFrame(const model::MCState& state, int frameNumber);
    void writeXYZFrame(const model::MCState& state, int frameNumber);
    void writeDATFrame(const model::MCState& state, int frameNumber);
    void writeTOPData(const model::MCState& state);
    
    std::string formatPDBAtom(int serial, const std::string& name, 
                              const std::string& resName, int resSeq,
                              double x, double y, double z,
                              double occupancy = 1.0, double tempFactor = 0.0,
                              const std::string& element = "");
    
    std::ofstream file_;
    Format format_;
    std::string filename_;
    int frameCount_;
    bool isOpen_;
};

// Utility class for writing analysis data
class DataWriter {
public:
    DataWriter(const std::string& filename);
    ~DataWriter();
    
    // Write header with column names
    void writeHeader(const std::vector<std::string>& columns);
    
    // Write data row
    void writeRow(const std::vector<double>& values);
    
    // Write comment
    void writeComment(const std::string& comment);
    
    void close();
    
private:
    std::ofstream file_;
    bool headerWritten_;
};

} // namespace io
} // namespace pygcmc