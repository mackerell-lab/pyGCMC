#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_HPP

#include "../../movement/gcmc/GCMCEngine.hpp"
#include "../../movement/gcmc/GCMCAcceptance.hpp"
#include "../../movement/gcmc/GCMCStatistics.hpp"
#include "../../movement/reservoir/fragment_reservoir.hpp"
#include "../../movement/reservoir/MultiTypeReservoir.hpp"
#include "../../energy/EnergyModule.hpp"
#include "../../../../model/param/ParamMain.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../../../io/parameters/InpParserMain.hpp"
#include "../../../../io/parameters/InpParserGCMC.hpp"
#include "../../../../system/log/LogMain.hpp"
#include "../stats/StatisticsTracker.hpp"
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

/**
 * @brief Main GCMC simulation class that coordinates all CPU platform components
 * 
 * This class provides a complete GCMC simulation framework that:
 * - Reads input files (INP format compatible with gcmc_gpu/gcmc_opencl)
 * - Manages multiple fragment types with individual chemical potentials
 * - Performs MC moves using GCMCEngine
 * - Tracks statistics and convergence
 * - Outputs trajectories and analysis
 */
class GCMCSimulation {
public:
    /**
     * @brief Configuration for the simulation
     */
    struct Config {
        std::string inputFile;              // Path to INP file
        std::string outputPrefix = "gcmc";  // Output file prefix
        int printFrequency = 1000;          // Statistics print frequency
        int trajectoryFrequency = 10000;    // Trajectory save frequency
        int checkpointFrequency = 100000;   // Checkpoint save frequency
        bool verbose = false;               // Verbose output
        int randomSeed = -1;                // Random seed (-1 for auto)
        
        // Performance options
        bool enableStatistics = true;       // Enable statistics collection
        int statisticsInterval = 1000;      // Statistics sampling interval
        bool storeProbabilities = false;    // Store acceptance probabilities
        
        // Advanced options
        bool enableAdaptiveSampling = false;  // Adjust move probabilities
        bool enableEnergyMinimization = false; // Minimize after insertion
        double convergenceTolerance = 0.01;    // Convergence criterion
    };
    
    /**
     * @brief Statistics for the simulation
     */
    struct Statistics {
        // Overall statistics
        int totalSteps = 0;
        int acceptedMoves = 0;
        double acceptanceRate = 0.0;
        
        // Move-specific statistics
        std::map<std::string, int> moveAttempts;
        std::map<std::string, int> moveAccepted;
        std::map<std::string, double> moveAcceptanceRates;
        
        // Fragment-specific statistics
        std::map<std::string, int> fragmentCounts;
        std::map<std::string, double> fragmentDensities;
        std::map<std::string, double> fragmentAcceptanceRates;
        
        // Energy statistics
        double currentEnergy = 0.0;
        double averageEnergy = 0.0;
        double energyStdDev = 0.0;
        std::vector<double> energyHistory;
        
        // Timing
        double totalTime = 0.0;  // seconds
        double timePerStep = 0.0;  // seconds
        double stepsPerSecond = 0.0;
    };
    
    /**
     * @brief Fragment type information
     */
    struct FragmentInfo {
        std::string name;
        int typeId;
        double concentration;      // Target concentration (M)
        double chemicalPotential;  // Chemical potential (kJ/mol)
        double activity;           // Activity (computed from μ)
        double probability;        // Selection probability
        int maxCount;             // Maximum number allowed
        int currentCount = 0;     // Current number in system
        
        // Template information
        movement::FragmentTemplate template_;
        
        // Statistics
        int insertAttempts = 0;
        int insertAccepted = 0;
        int deleteAttempts = 0;
        int deleteAccepted = 0;
    };
    
    // Constructor and destructor
    explicit GCMCSimulation(const Config& config);
    ~GCMCSimulation();
    
    // Main simulation methods
    bool initialize();           // Load input and setup system
    bool run();                 // Run the simulation
    void finalize();            // Clean up and write final output
    
    // Control methods
    void stop() { running_ = false; }
    bool isRunning() const { return running_; }
    
    // Analysis methods
    Statistics getStatistics() const { return stats_; }
    void printStatistics() const;
    void saveTrajectory(const std::string& filename) const;
    void saveCheckpoint(const std::string& filename) const;
    bool loadCheckpoint(const std::string& filename);
    
    // Configuration access
    Config getConfig() const { return config_; }
    void updateConfig(const Config& config) { config_ = config; }
    
    // Fragment information access
    std::vector<FragmentInfo> getFragmentInfo() const { return fragmentTypes_; }
    
private:
    // Configuration
    Config config_;
    bool initialized_ = false;
    bool running_ = false;
    
    // Core components
    std::unique_ptr<model::param::Param> params_;
    std::unique_ptr<model::montecarlo::MCState> state_;
    std::unique_ptr<movement::gcmc::GCMCEngine> engine_;
    std::unique_ptr<movement::gcmc::GCMCAcceptance> acceptance_;
    std::unique_ptr<movement::MultiTypeReservoir> reservoir_;
    movement::gcmc::GCMCStatistics statistics_;
    
    // Fragment management
    std::vector<FragmentInfo> fragmentTypes_;
    std::map<std::string, int> fragmentNameToId_;
    
    // Statistics
    Statistics stats_;
    StatisticsTracker simulationStats_;  // New modular statistics tracker
    std::chrono::steady_clock::time_point startTime_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;
    
    // Helper methods
    bool loadParameters();
    bool setupSystem();
    bool setupFragments();
    bool setupAcceptance();
    bool setupEngine();
    bool performMCStep();
    
    // Move selection
    enum MoveType {
        INSERT,
        DELETE,
        TRANSLATE,
        ROTATE
    };
    MoveType selectMoveType();
    int selectFragmentType();
    int selectActiveFragment();
    
    // Energy and analysis
    double calculateSystemEnergy();
    void updateStatistics();
    bool checkConvergence();
    
    // Output methods
    void writeStatistics(int step);
    void writeTrajectory(int step);
    void writeCheckpoint(int step);
    void writeFinalResults();
    
    // Logging helper
    template<typename... Args>
    void log(const std::string& format, Args... args) const;
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_HPP