#pragma once

/**
 * @file SimulationAPI.hpp
 * @brief Public API for GCMC simulation management
 *
 * This file provides the unified public interface for GCMC simulations,
 * consolidating functionality from various simulation modules.
 */

#include "../../../model/ModelModule.hpp"
#include "../../../system/common/SystemLogger.hpp"
#include "core/SimulationCore.hpp"
#include "setup/SystemInitializer.hpp"
#include "stats/StatisticsTracker.hpp"
#include "io/SimulationIO.hpp"
#include "gcmc/GCMCController.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

using namespace pygcmc::system::common;
using namespace pygcmc::model::montecarlo;

// Re-export key classes for convenience
using SimulationCore = core::SimulationCore;
using SystemInitializer = setup::SystemInitializer;
using StatisticsTracker = stats::StatisticsTracker;
using SimulationIO = io::SimulationIO;
using GCMCController = gcmc::GCMCController;

/**
 * @brief Main GCMC Simulation interface
 *
 * This class provides a unified interface for GCMC simulations,
 * delegating implementation to specialized modules.
 */
class GCMCSimulation {
public:
    /**
     * @brief Configuration for the simulation
     */
    struct Config {
        // Input/Output configuration
        std::string inputFile;              // Path to INP file
        std::string outputPrefix = "gcmc";  // Output file prefix

        // Simulation parameters
        int numSteps = 100000;               // Total number of MC steps
        int equilibrationSteps = 10000;      // Equilibration steps
        double temperature = 298.15;          // Temperature in K

        // Output frequencies
        int printFrequency = 1000;           // Statistics print frequency
        int trajectoryFrequency = 10000;     // Trajectory save frequency
        int checkpointFrequency = 100000;    // Checkpoint save frequency

        // Control options
        bool verbose = false;                // Verbose output
        int randomSeed = -1;                 // Random seed (-1 for auto)

        // Performance options
        bool enableStatistics = true;        // Enable statistics collection
        int statisticsInterval = 1000;       // Statistics sampling interval

        // Advanced options
        bool enableAdaptiveSampling = false; // Adjust move probabilities
        double convergenceTolerance = 0.01;  // Convergence criterion

        // Energy method
        std::string energyMethod = "direct"; // direct, ewald, pme, pgp
        double cutoff = 12.0;                // Cutoff distance (Angstrom)

        // Move probabilities (should sum to 1.0)
        double insertionProb = 0.25;
        double deletionProb = 0.25;
        double translationProb = 0.25;
        double rotationProb = 0.25;
    };

    // Constructor and destructor
    explicit GCMCSimulation(const Config& config = Config());
    ~GCMCSimulation();

    // Prevent copying
    GCMCSimulation(const GCMCSimulation&) = delete;
    GCMCSimulation& operator=(const GCMCSimulation&) = delete;

    // Allow moving
    GCMCSimulation(GCMCSimulation&&) = default;
    GCMCSimulation& operator=(GCMCSimulation&&) = default;

    // ========================================================================
    // Main simulation control
    // ========================================================================

    /**
     * @brief Initialize simulation from configuration
     * @return true if successful
     */
    bool initialize();

    /**
     * @brief Initialize simulation from INP file
     * @param inpFile Path to INP file
     * @return true if successful
     */
    bool initializeFromFile(const std::string& inpFile);

    /**
     * @brief Run the simulation
     * @param nSteps Number of steps to run (0 = use config.numSteps)
     * @return true if successful
     */
    bool run(int nSteps = 0);

    /**
     * @brief Finalize simulation and write output
     */
    void finalize();

    /**
     * @brief Stop running simulation
     */
    void stop();

    /**
     * @brief Check if simulation is running
     */
    bool isRunning() const;

    // ========================================================================
    // State access
    // ========================================================================

    /**
     * @brief Get the current MC state
     */
    MCState& getState();
    const MCState& getState() const;

    /**
     * @brief Get simulation configuration
     */
    const Config& getConfig() const;
    void updateConfig(const Config& config);

    // ========================================================================
    // Statistics and analysis
    // ========================================================================

    /**
     * @brief Get current statistics
     */
    const StatisticsTracker& getStatistics() const;

    /**
     * @brief Print current statistics
     */
    void printStatistics() const;

    /**
     * @brief Export statistics to file
     */
    void exportStatistics(const std::string& filename) const;

    /**
     * @brief Get energy components
     */
    std::pair<double, double> getEnergyComponents() const;

    /**
     * @brief Get molecule counts by type
     */
    std::map<std::string, int> getMoleculeCounts() const;

    // ========================================================================
    // Trajectory and checkpoint management
    // ========================================================================

    /**
     * @brief Write current frame to trajectory
     */
    void writeTrajectory(int step = -1);

    /**
     * @brief Write checkpoint file
     */
    void writeCheckpoint(const std::string& filename = "");

    /**
     * @brief Load checkpoint file
     */
    bool loadCheckpoint(const std::string& filename);

    // ========================================================================
    // Fragment management
    // ========================================================================

    /**
     * @brief Add fragment type to simulation
     */
    void addFragmentType(const std::string& name,
                        const std::vector<MCAtom>& atoms,
                        double chemicalPotential = -15.7,
                        double targetDensity = 0.0);

    /**
     * @brief Set fragment reservoir from file
     */
    bool setFragmentReservoir(const std::string& fragmentFile);

    /**
     * @brief Get fragment information
     */
    std::vector<SystemInitializer::FragmentConfig> getFragmentInfo() const;

    // ========================================================================
    // Move control
    // ========================================================================

    /**
     * @brief Perform a single MC move
     */
    bool performMove();

    /**
     * @brief Set move probabilities
     */
    void setMoveProbabilities(double insertion, double deletion,
                             double translation, double rotation);

    /**
     * @brief Get current move probabilities
     */
    std::map<std::string, double> getMoveProbabilities() const;

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

// ============================================================================
// Convenience Functions
// ============================================================================

/**
 * @brief Create and run a simple GCMC simulation
 */
inline void runSimpleGCMC(const std::string& inpFile, int nSteps) {
    if (SystemLogger::isInfoEnabled()) {
        SystemLogger::info("Running simple GCMC simulation from ", inpFile);
    }

    GCMCSimulation::Config config;
    config.inputFile = inpFile;
    config.numSteps = nSteps;

    GCMCSimulation sim(config);
    if (sim.initialize()) {
        sim.run();
        sim.finalize();
    }
}

/**
 * @brief Create a GCMC simulation with custom configuration
 */
inline std::unique_ptr<GCMCSimulation> createGCMCSimulation(
    const GCMCSimulation::Config& config) {
    return std::make_unique<GCMCSimulation>(config);
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // Include guard
