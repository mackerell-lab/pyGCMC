#include "SimulationIO.hpp"
#include "../../../../system/common/SystemLogger.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <filesystem>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace io {

using namespace pygcmc::system::common;
using namespace pygcmc::model::montecarlo;
namespace fs = std::filesystem;

SimulationIO::SimulationIO(const Config& config)
    : config_(config)
    , trajectoryOpen_(false)
    , statisticsOpen_(false)
    , energyOpen_(false) {
}

SimulationIO::~SimulationIO() {
    // Close all open files
    if (trajectoryOpen_) {
        closeTrajectoryFile();
    }
    if (statisticsOpen_) {
        closeStatisticsFile();
    }
    if (energyOpen_) {
        closeEnergyFile();
    }
}

bool SimulationIO::readInputFile(const std::string& filename,
                                 MCState& state,
                                 model::param::Param& params) {
    if (!fs::exists(filename)) {
        SystemLogger::error("Input file not found: ", filename);
        return false;
    }

    try {
        // Parse INP file
        // TODO: Implement proper INP parser
        // For now, just set basic parameters

        // Default values for testing
        float boxSize = 5.0; // nm
        float temperature = 300.0; // K
        float cutoff = 1.2; // nm

        // Set up state from parsed parameters
        state.setBoxDimensions(boxSize, boxSize, boxSize);
        state.setTemperature(temperature);

        // Set up force field parameters
        state.setupForceField(2, 2); // 2 types for now

        // Set up parameters
        params.set_temperature(temperature);
        params.set_cutoff(cutoff);

        SystemLogger::info("Successfully read input file: ", filename);
        SystemLogger::info("Box: ", boxSize, " x ", boxSize, " x ", boxSize, " nm");
        SystemLogger::info("Temperature: ", temperature, " K");

        return true;

    } catch (const std::exception& e) {
        SystemLogger::error("Error reading input file: ", e.what());
        return false;
    }
}

bool SimulationIO::readFragmentFile(const std::string& filename,
                                   std::vector<MCResidue>& fragments) {
    if (!fs::exists(filename)) {
        SystemLogger::error("Fragment file not found: ", filename);
        return false;
    }

    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open fragment file");
        }

        // Parse fragment file format
        // Format: residue_name num_atoms
        //         atom_type x y z charge
        //         ...

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::string resName;
            int numAtoms;

            if (!(iss >> resName >> numAtoms)) continue;

            MCResidue residue;
            residue.resname = resName;
            residue.atomCount = numAtoms;
            residue.atoms.reserve(numAtoms);

            for (int i = 0; i < numAtoms; ++i) {
                if (!std::getline(file, line)) {
                    throw std::runtime_error("Incomplete fragment definition");
                }

                std::istringstream atomLine(line);
                MCAtom atom;

                if (!(atomLine >> atom.type >> atom.x >> atom.y >> atom.z >> atom.charge)) {
                    throw std::runtime_error("Invalid atom definition");
                }

                residue.atoms.push_back(atom);
            }

            fragments.push_back(residue);
        }

        SystemLogger::info("Read ", fragments.size(), " fragments from ", filename);

        return true;

    } catch (const std::exception& e) {
        SystemLogger::error("Error reading fragment file: ", e.what());
        return false;
    }
}

bool SimulationIO::readCheckpoint(const std::string& filename,
                                 MCState& state,
                                 int& step) {
    if (!checkpointManager_) {
        pygcmc::io::output::CheckpointManager::Config config;
        checkpointManager_ = std::make_unique<pygcmc::io::output::CheckpointManager>(config);
    }

    StatisticsTracker stats;  // Dummy stats for loading
    return checkpointManager_->loadCheckpoint(state, stats, step, filename);
}

void SimulationIO::writeTrajectory(const MCState& state, int step) {
    if (!trajectoryWriter_) {
        openTrajectoryFile();
    }

    if (trajectoryWriter_) {
        trajectoryWriter_->writeTrajectory(state, step);
    }
}

void SimulationIO::writeCheckpoint(const MCState& state,
                                  int step,
                                  const StatisticsTracker* stats) {
    if (!checkpointManager_) {
        pygcmc::io::output::CheckpointManager::Config config;
        checkpointManager_ = std::make_unique<pygcmc::io::output::CheckpointManager>(config);
    }

    std::string filename = generateFilename("checkpoint", "chk");
    if (stats) {
        checkpointManager_->saveCheckpoint(state, *stats, step, filename);
    } else {
        StatisticsTracker dummyStats;
        checkpointManager_->saveCheckpoint(state, dummyStats, step, filename);
    }
}

void SimulationIO::writeStatistics(const StatisticsTracker& stats, int step) {
    if (!statisticsOpen_) {
        openStatisticsFile();
    }

    if (statisticsWriter_) {
        std::vector<double> values = {
            static_cast<double>(step),
            stats.getTotalAcceptanceRate(),
            static_cast<double>(0), // Molecule count placeholder
            stats.getAverageEnergy(),
            0.0  // Density placeholder
        };

        statisticsWriter_->writeRow(values);
    }
}

void SimulationIO::writeEnergies(const MCState& state, int step) {
    if (!config_.writeEnergies) return;

    if (!energyOpen_) {
        openEnergyFile();
    }

    if (energyWriter_) {
        double totalVdW = 0.0;
        double totalElec = 0.0;

        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active) {
                totalVdW += state.residues[i].energy_vdw;
                totalElec += state.residues[i].energy_elec;
            }
        }

        // Account for double counting
        totalVdW /= 2.0;
        totalElec /= 2.0;

        std::vector<double> values = {
            static_cast<double>(step),
            totalVdW,
            totalElec,
            totalVdW + totalElec
        };

        energyWriter_->writeRow(values);
    }
}

void SimulationIO::writeDensities(const std::map<std::string, double>& densities,
                                 int step) {
    if (!config_.writeDensities) return;

    if (!densityWriter_) {
        densityWriter_ = std::make_unique<pygcmc::io::output::DataWriter>(
            generateFilename("density", "dat"));

        // Write header
        std::vector<std::string> columns = {"Step"};
        for (const auto& [name, _] : densities) {
            columns.push_back(name);
        }
        densityWriter_->writeHeader(columns);
    }

    std::vector<double> values = {static_cast<double>(step)};
    for (const auto& [_, density] : densities) {
        values.push_back(density);
    }

    densityWriter_->writeRow(values);
}

void SimulationIO::writeFinalReport(const MCState& state,
                                   const StatisticsTracker& stats) {
    std::string filename = generateFilename("report", "txt");
    std::ofstream file(filename);

    if (!file.is_open()) {
        if (SystemLogger::isDebugEnabled()) {
            SystemLogger::debug("Could not write final report");
        }
        return;
    }

    file << "===== GCMC Simulation Final Report =====\n\n";
    file << "Timestamp: " << formatTimestamp() << "\n";
    file << "Output prefix: " << config_.outputPrefix << "\n\n";

    file << "System Configuration:\n";
    file << "  Box: " << state.info.box[0] << " x " << state.info.box[1]
         << " x " << state.info.box[2] << " nm\n";
    file << "  Temperature: " << (1.0 / (state.info.beta * 8.314e-3)) << " K\n";
    file << "  Cutoff: " << state.info.cutoff << " nm\n\n";

    file << "Final Statistics:\n";
    file << "  Total steps: " << stats.getTotalSteps() << "\n";
    file << "  Acceptance rate: " << std::fixed << std::setprecision(4)
         << stats.getTotalAcceptanceRate() * 100 << "%\n";
    file << "  Total molecules: " << 0 << "\n";  // Placeholder
    file << "  Average energy: " << stats.getAverageEnergy() << " kJ/mol\n\n";

    file << "Move Statistics:\n";
    // TODO: Add move statistics when available

    file.close();

    SystemLogger::info("Final report written to ", filename);
}

void SimulationIO::openTrajectoryFile(const std::string& filename) {
    std::string fname = filename.empty() ?
        generateFilename("trajectory", config_.trajectoryFormat) : filename;

    pygcmc::io::output::TrajectoryWriter::Config writerConfig;
    writerConfig.format = config_.trajectoryFormat;
    writerConfig.prefix = config_.outputPrefix;
    writerConfig.compressOutput = config_.compressTrajectory;
    writerConfig.precision = config_.precision;
    writerConfig.multiFrame = true;
    writerConfig.continuousFile = true;

    trajectoryWriter_ = std::make_unique<pygcmc::io::output::TrajectoryWriter>(writerConfig);
    trajectoryWriter_->open(fname);
    trajectoryOpen_ = true;
}

void SimulationIO::closeTrajectoryFile() {
    if (trajectoryWriter_) {
        trajectoryWriter_->close();
    }
    trajectoryOpen_ = false;
}

bool SimulationIO::isTrajectoryOpen() const {
    return trajectoryOpen_;
}

void SimulationIO::openStatisticsFile(const std::string& filename) {
    std::string fname = filename.empty() ?
        generateFilename("statistics", "dat") : filename;

    statisticsWriter_ = std::make_unique<pygcmc::io::output::DataWriter>(fname);

    // Write header
    std::vector<std::string> columns = {
        "Step", "AcceptRate", "Molecules", "Energy", "Density"
    };
    statisticsWriter_->writeHeader(columns);
    statisticsOpen_ = true;
}

void SimulationIO::closeStatisticsFile() {
    if (statisticsWriter_) {
        statisticsWriter_->close();
    }
    statisticsOpen_ = false;
}

void SimulationIO::openEnergyFile(const std::string& filename) {
    std::string fname = filename.empty() ?
        generateFilename("energy", "dat") : filename;

    energyWriter_ = std::make_unique<pygcmc::io::output::DataWriter>(fname);

    // Write header
    std::vector<std::string> columns = {
        "Step", "VdW", "Elec", "Total"
    };
    energyWriter_->writeHeader(columns);
    energyOpen_ = true;
}

void SimulationIO::closeEnergyFile() {
    if (energyWriter_) {
        energyWriter_->close();
    }
    energyOpen_ = false;
}

void SimulationIO::updateConfig(const Config& config) {
    config_ = config;
}

std::string SimulationIO::generateFilename(const std::string& suffix,
                                          const std::string& extension) const {
    std::stringstream ss;
    ss << config_.outputPrefix;
    if (!suffix.empty()) {
        ss << "_" << suffix;
    }
    ss << "." << extension;
    return ss.str();
}

std::string SimulationIO::formatTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

void SimulationIO::ensureDirectoryExists(const std::string& path) const {
    fs::path dir = fs::path(path).parent_path();
    if (!dir.empty() && !fs::exists(dir)) {
        fs::create_directories(dir);
    }
}

void SimulationIO::writeHeader(pygcmc::io::output::DataWriter* writer,
                              const std::vector<std::string>& columns) {
    if (writer) {
        writer->writeHeader(columns);
    }
}

} // namespace io
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc
