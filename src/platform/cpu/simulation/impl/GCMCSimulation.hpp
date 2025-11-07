#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_HPP

#include "../../movement/gcmc/GCMCEngine.hpp"
#include "../../movement/gcmc/GCMCAcceptance.hpp"
#include "../../movement/gcmc/GCMCStatistics.hpp"
#include "../../movement/reservoir/fragment_reservoir.hpp"
#include "../../movement/reservoir/MultiTypeReservoir.hpp"
#include "../../movement/bias/CavityBias.hpp"
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
#include <array>

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
        int printFrequency = 0;             // Statistics print frequency (0 = use INP nprint)
        int trajectoryFrequency = 10000;    // Trajectory save frequency
        int checkpointFrequency = 0;         // Checkpoint save frequency (0 = disabled)
        bool verbose = false;               // Verbose output
        int randomSeed = -1;                // Random seed (-1 for auto)
        
        // Performance options
        bool enableStatistics = true;       // Enable statistics collection
        int movesPerStep = 1;               // Number of moves per MC step (default 1 for proper GCMC)
        int statisticsInterval = 1000;      // Statistics sampling interval
        bool storeProbabilities = false;    // Store acceptance probabilities
        
        // Advanced options
        bool enableAdaptiveSampling = false;  // Adjust move probabilities
        bool enableEnergyMinimization = false; // Minimize after insertion
        double convergenceTolerance = 0.01;    // Convergence criterion
        int maxMoleculesPerType = 10000;       // Max molecules per fragment type (-1 = disabled)
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
     * @brief Detailed acceptance counters for diagnostics
     */
    struct Counters {
        uint64_t attemptsInsert = 0;
        uint64_t acceptsInsert = 0;
        uint64_t attemptsDelete = 0;
        uint64_t acceptsDelete = 0;
        uint64_t attemptsTranslate = 0;
        uint64_t acceptsTranslate = 0;
        uint64_t attemptsRotate = 0;
        uint64_t acceptsRotate = 0;

        // Since last print
        uint64_t attemptsSinceLastPrint = 0;
        uint64_t acceptsSinceLastPrint = 0;

        double getInsDelOverallRate() const {
            uint64_t total = attemptsInsert + attemptsDelete;
            return total > 0 ?
                static_cast<double>(acceptsInsert + acceptsDelete) / total : 0.0;
        }
    };

    /**
     * @brief Detailed acceptance diagnostics record
     */
    struct AcceptanceRecord {
        enum MoveType { INSERT, DELETE, TRANSLATE, ROTATE };

        MoveType moveType;
        int species;
        int nBefore;       // Number of molecules of this species before the move
        int cbmcTrials;    // Number of CBMC trial configurations (K)
        int step;
        double deltaU;
        double betaDeltaU;
        double mu;
        double betaMu;
        double z;  // Activity
        double qForward;   // Forward proposal probability
        double qReverse;   // Reverse proposal probability
        double proposalRatio;
        double vEff;       // Effective volume
        double vBox = 0.0; // Base box volume at time of move
        double cavityFraction = 1.0; // Recorded cavity fraction (f_cav)
        double rosenbluthWeight = 1.0; // Recorded Rosenbluth ratio (W/K or K/W)
        double pAcc;       // Calculated acceptance probability
        double bias = 1.0; // Total bias factor applied in acceptance
        double u;          // Random number used
        bool accepted;

        // For biased moves
        double wForward = 1.0;
        double wReverse = 1.0;
        double wCavity = 1.0;   // Cavity bias factor (V_cavity / V_box)
    };

    /**
     * @brief Fragment type information
     */
    struct FragmentInfo {
        std::string name;
        int typeId = 0;
        double concentration = 0.0;      // Target concentration (M)
        double chemicalPotential = 0.0;  // Chemical potential (kJ/mol)
        double activity = 0.0;           // Activity (computed from μ)
        double probability = 0.0;        // Selection probability
        int maxCount = 0;                // Maximum number allowed
        int currentCount = 0;     // Current number in system
        int confBiasTrials = 1;   // Number of configuration bias trials

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
    void saveTopology(const std::string& filename) const;
    void saveCheckpoint(const std::string& filename) const;
    bool loadCheckpoint(const std::string& filename);
    
    // Configuration access
    Config getConfig() const { return config_; }
    void updateConfig(const Config& config) { config_ = config; }
    
    // Fragment information access
    std::vector<FragmentInfo> getFragmentInfo() const { return fragmentTypes_; }

    // Diagnostic methods
    void enableDiagnostics(size_t bufferSize = 4096);
    bool isDiagnosticsEnabled() const { return diagnosticsEnabled_; }
    Counters getCounters() const { return counters_; }
    AcceptanceRecord getLastMove() const;
    std::vector<AcceptanceRecord> getMoves(size_t n) const;
    void dumpLJMatrix() const;
    void dumpAcceptanceLog(const std::string& filename) const;

private:
    // Statistics output
    void writeStatisticsDAT(int step);
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
    std::unique_ptr<movement::CavityManager> cavityManager_;  // Cavity bias manager
    movement::RegionConstraint* regionConstraint_ = nullptr;  // Raw pointer to region constraint (owned by engine)
    movement::gcmc::GCMCStatistics statistics_;
    
    // Fragment management
    std::vector<FragmentInfo> fragmentTypes_;
    std::map<std::string, int> fragmentNameToId_;
    std::map<std::string, movement::FragmentTemplate> fragmentTemplatesFromBuilder_;
    std::map<std::string, size_t> atomTypeNameToIndex_;  // Atom type name to force field index mapping
    
    // Force field from builder (if loaded)
    std::shared_ptr<model::ForceField> forceFieldFromBuilder_;
    
    // Statistics
    Statistics stats_;
    StatisticsTracker simulationStats_;  // New modular statistics tracker
    std::chrono::steady_clock::time_point startTime_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;

    // Proposal probabilities for detailed balance with target_numwaters
    double lastProposalPInsert_ = 0.25;
    double lastProposalPDelete_ = 0.25;
    double currentProposalRatio_ = 1.0;

    // Diagnostics
    bool diagnosticsEnabled_ = false;
    Counters counters_;
    std::vector<AcceptanceRecord> acceptanceBuffer_;
    size_t bufferSize_ = 4096;
    size_t bufferIndex_ = 0;

    // Helper methods
    bool loadParameters();
    void printParameterSummary();
    bool setupSystem();
    bool setupFragments();
    bool setupAcceptance();
    bool setupEngine();
    bool performMCStep();
    bool performSingleMove();  // Performs a single GCMC move

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
    void outputWaterDensity(int step);
    
    // Logging helper
    template<typename... Args>
    void log(const std::string& format, Args... args) const;

    // === Added: per-fragment move probabilities (Ins/Del/Trn/Rot) ===
    // CDF per fragment: [P(Ins), P(Del), P(Trn), P(Rot)] cumulative
    std::vector<std::array<double, 4>> fragmentMoveCDF_;
    void buildPerFragmentMoveCDF();
    MoveType selectMoveForFragment(int fragType);

    // === Added: nbar modes integration ===
    void updateActivitiesForNbar();
    static bool isWaterName(const std::string& name);
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_GCMC_SIMULATION_HPP
