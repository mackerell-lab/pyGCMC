#include "CheckpointManager.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

CheckpointManager::CheckpointManager(const Config& config)
    : config_(config) {
}

CheckpointManager::~CheckpointManager() {
}

bool CheckpointManager::saveCheckpoint(const model::montecarlo::MCState& state,
                                      const StatisticsTracker& stats,
                                      int step,
                                      const std::string& filename) {
    std::string fname = filename.empty() ? generateFilename(step) : filename;
    
    std::ofstream out(fname, std::ios::binary);
    if (!out) {
        std::cerr << "Failed to open checkpoint file: " << fname << std::endl;
        return false;
    }
    
    // Write header
    if (!writeHeader(out, config_.version)) {
        return false;
    }
    
    // Write step
    out.write(reinterpret_cast<const char*>(&step), sizeof(step));
    
    // Write state dimensions
    out.write(reinterpret_cast<const char*>(&state.activeAtomCount), 
             sizeof(state.activeAtomCount));
    out.write(reinterpret_cast<const char*>(&state.activeResidueCount), 
             sizeof(state.activeResidueCount));
    
    // Write box
    for (double dim : state.periodicBox) {
        out.write(reinterpret_cast<const char*>(&dim), sizeof(dim));
    }
    
    // Write atoms (simplified)
    for (int i = 0; i < state.activeAtomCount; ++i) {
        const auto& atom = state.atoms[i];
        out.write(reinterpret_cast<const char*>(&atom.x), sizeof(atom.x));
        out.write(reinterpret_cast<const char*>(&atom.y), sizeof(atom.y));
        out.write(reinterpret_cast<const char*>(&atom.z), sizeof(atom.z));
    }
    
    // Write statistics (basic info)
    size_t totalSteps = stats.getTotalSteps();
    size_t totalAccepted = stats.getTotalAccepted();
    double currentEnergy = stats.getCurrentEnergy();
    
    out.write(reinterpret_cast<const char*>(&totalSteps), sizeof(totalSteps));
    out.write(reinterpret_cast<const char*>(&totalAccepted), sizeof(totalAccepted));
    out.write(reinterpret_cast<const char*>(&currentEnergy), sizeof(currentEnergy));
    
    return true;
}

bool CheckpointManager::loadCheckpoint(model::montecarlo::MCState& state,
                                      StatisticsTracker& stats,
                                      int& step,
                                      const std::string& filename) {
    if (!validateCheckpoint(filename)) {
        return false;
    }
    
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        std::cerr << "Failed to open checkpoint file: " << filename << std::endl;
        return false;
    }
    
    // Read header
    int version;
    if (!readHeader(in, version)) {
        return false;
    }
    
    // Read step
    in.read(reinterpret_cast<char*>(&step), sizeof(step));
    
    // Read state dimensions
    in.read(reinterpret_cast<char*>(&state.activeAtomCount), 
           sizeof(state.activeAtomCount));
    in.read(reinterpret_cast<char*>(&state.activeResidueCount), 
           sizeof(state.activeResidueCount));
    
    // Read box
    state.periodicBox.resize(3);
    for (double& dim : state.periodicBox) {
        in.read(reinterpret_cast<char*>(&dim), sizeof(dim));
    }
    
    // Read atoms
    state.atoms.resize(state.activeAtomCount);
    for (int i = 0; i < state.activeAtomCount; ++i) {
        auto& atom = state.atoms[i];
        in.read(reinterpret_cast<char*>(&atom.x), sizeof(atom.x));
        in.read(reinterpret_cast<char*>(&atom.y), sizeof(atom.y));
        in.read(reinterpret_cast<char*>(&atom.z), sizeof(atom.z));
    }
    
    // Read statistics
    size_t totalSteps, totalAccepted;
    double currentEnergy;
    
    in.read(reinterpret_cast<char*>(&totalSteps), sizeof(totalSteps));
    in.read(reinterpret_cast<char*>(&totalAccepted), sizeof(totalAccepted));
    in.read(reinterpret_cast<char*>(&currentEnergy), sizeof(currentEnergy));
    
    // Reset stats and restore basic info
    stats.reset();
    for (size_t i = 0; i < totalSteps; ++i) {
        stats.recordMove("restore", "fragment", i < totalAccepted);
    }
    stats.recordEnergy(currentEnergy);
    
    return true;
}

bool CheckpointManager::validateCheckpoint(const std::string& filename) const {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        return false;
    }
    
    // Check magic number
    char magic[5] = {0};
    in.read(magic, 4);
    
    return std::string(magic) == "GCMC";
}

std::string CheckpointManager::generateFilename(int step) const {
    std::stringstream ss;
    ss << config_.prefix << "_";
    if (step >= 0) {
        ss << std::setfill('0') << std::setw(8) << step << "_";
    }
    ss << "checkpoint.chk";
    return ss.str();
}

bool CheckpointManager::writeHeader(std::ofstream& out, int version) const {
    const char* magic = "GCMC";
    out.write(magic, 4);
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    return out.good();
}

bool CheckpointManager::readHeader(std::ifstream& in, int& version) const {
    char magic[5] = {0};
    in.read(magic, 4);
    if (std::string(magic) != "GCMC") {
        return false;
    }
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    return in.good();
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc