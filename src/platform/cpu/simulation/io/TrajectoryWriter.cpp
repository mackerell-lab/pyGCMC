#include "TrajectoryWriter.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

TrajectoryWriter::TrajectoryWriter(const Config& config)
    : config_(config) {
}

TrajectoryWriter::~TrajectoryWriter() {
}

bool TrajectoryWriter::writeTrajectory(const model::montecarlo::MCState& state,
                                      int step,
                                      const std::string& filename) {
    std::string fname = filename.empty() ? generateFilename(step) : filename;
    
    std::ofstream out(fname);
    if (!out) {
        std::cerr << "Failed to open trajectory file: " << fname << std::endl;
        return false;
    }
    
    if (config_.format == "pdb") {
        return writePDB(state, step, out);
    } else if (config_.format == "xyz") {
        return writeXYZ(state, step, out);
    }
    
    return false;
}

bool TrajectoryWriter::writePDB(const model::montecarlo::MCState& state,
                               int step,
                               std::ofstream& out) {
    writePDBHeader(out, step);
    writePDBBox(out, state.periodicBox);
    
    // Write atoms
    for (int i = 0; i < state.activeAtomCount; ++i) {
        const auto& atom = state.atoms[i];
        
        std::string resName = "UNK";
        int resSeq = i + 1;  // Simple sequential numbering
        
        out << formatPDBAtom(i + 1, atom.name, resName, resSeq,
                           atom.x * 10.0,  // nm to Angstrom
                           atom.y * 10.0,
                           atom.z * 10.0) << std::endl;
    }
    
    out << "END" << std::endl;
    return true;
}

bool TrajectoryWriter::writeXYZ(const model::montecarlo::MCState& state,
                               int step,
                               std::ofstream& out) {
    // Simple XYZ format
    out << state.activeAtomCount << std::endl;
    out << "Step " << step << std::endl;
    
    for (int i = 0; i < state.activeAtomCount; ++i) {
        const auto& atom = state.atoms[i];
        out << atom.name << " "
            << atom.x * 10.0 << " "
            << atom.y * 10.0 << " "
            << atom.z * 10.0 << std::endl;
    }
    
    return true;
}

std::string TrajectoryWriter::generateFilename(int step) const {
    std::stringstream ss;
    ss << config_.prefix << "_";
    ss << std::setfill('0') << std::setw(8) << step;
    
    if (config_.format == "pdb") {
        ss << ".pdb";
    } else if (config_.format == "xyz") {
        ss << ".xyz";
    } else {
        ss << ".dat";
    }
    
    return ss.str();
}

void TrajectoryWriter::writePDBHeader(std::ofstream& out, int step) const {
    out << "REMARK GCMC Trajectory" << std::endl;
    out << "REMARK Step: " << step << std::endl;
    
    std::time_t now = std::time(nullptr);
    out << "REMARK Generated: " 
        << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") 
        << std::endl;
}

void TrajectoryWriter::writePDBBox(std::ofstream& out,
                                  const std::vector<double>& box) const {
    if (box.size() >= 3) {
        out << "CRYST1" 
            << std::fixed << std::setprecision(3)
            << std::setw(9) << box[0] * 10.0  // nm to Angstrom
            << std::setw(9) << box[1] * 10.0
            << std::setw(9) << box[2] * 10.0
            << std::setw(7) << "90.00"
            << std::setw(7) << "90.00"
            << std::setw(7) << "90.00"
            << " P 1           1" << std::endl;
    }
}

std::string TrajectoryWriter::formatPDBAtom(int serial,
                                           const std::string& atomName,
                                           const std::string& resName,
                                           int resSeq,
                                           double x, double y, double z) const {
    std::stringstream ss;
    ss << std::left << std::setw(6) << "ATOM"
       << std::right << std::setw(5) << serial << " "
       << std::left << std::setw(4) << atomName
       << std::right << std::setw(4) << resName << " "
       << std::setw(4) << resSeq << "    "
       << std::fixed << std::setprecision(config_.precision)
       << std::setw(8) << x
       << std::setw(8) << y
       << std::setw(8) << z
       << std::setw(6) << std::setprecision(2) << 1.00
       << std::setw(6) << 0.00;
    
    return ss.str();
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc