#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_MODULE_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_MODULE_HPP

/**
 * @file GCMCModule.hpp
 * @brief GCMC Module - Complete Grand Canonical Monte Carlo Implementation
 * 
 * This module provides a comprehensive GCMC implementation that combines:
 * - Fragment Reservoir management
 * - Cavity-biased insertion
 * - Configurational bias (CBMC)
 * - Move selection and execution
 * - Statistical analysis
 * 
 * All functionality is implemented in C++ for maximum performance.
 * Python bindings are only for testing and validation.
 */

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../energy/common/EnergyInterface.hpp"
#include <memory>
#include <map>
#include <string>
#include <vector>

// Forward declarations to avoid circular dependencies
namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
    class FragmentReservoir;
    class FragmentTemplate;
    class CavityManager;
    class ConfigBiasManager;
namespace gcmc {
    class GCMCEngine;
    class GCMCMoveSelector;
    class GCMCBias;
    class GCMCAcceptance;
    class GCMCStats;
} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Main GCMC module class that orchestrates all components
 * 
 * This class provides a high-level interface for GCMC simulations,
 * managing all subcomponents and coordinating their interactions.
 */
class GCMCModule {
public:
    // Configuration structure
    struct Config {
        // Temperature and thermodynamics
        double temperature;
        double pressure;
        
        // Simulation parameters
        int equilibrationSteps;
        int productionSteps;
        int saveFrequency;
        
        // Move probabilities (should sum to 1.0)
        double insertProb;
        double deleteProb;
        double translateProb;
        double rotateProb;
        double swapProb;          // For multi-component
        
        // Bias parameters
        bool useCavityBias;
        double gridSpacing;        // nm
        double probeRadius;       // nm
        bool useConfigBias;
        int configTrials;
        
        // Advanced options
        bool useRegrowth;        // Regrowth moves
        bool useClusterMoves;    // Cluster translation/rotation
        double clusterCutoff;     // nm
        
        // Energy calculation
        platform::cpu::EnergyMethod energyMethod;
        double cutoff;             // nm
        
        // Output
        bool verbose;
        std::string trajectoryFile;
        std::string statisticsFile;
        
        // Constructor with default values
        Config()
            : temperature(300.0),
              pressure(1.0),
              equilibrationSteps(10000),
              productionSteps(100000),
              saveFrequency(1000),
              insertProb(0.2),
              deleteProb(0.2),
              translateProb(0.2),
              rotateProb(0.2),
              swapProb(0.2),
              useCavityBias(true),
              gridSpacing(0.2),
              probeRadius(0.14),
              useConfigBias(true),
              configTrials(10),
              useRegrowth(false),
              useClusterMoves(false),
              clusterCutoff(0.35),
              energyMethod(platform::cpu::EnergyMethod::PME),
              cutoff(1.2),
              verbose(false),
              trajectoryFile(""),
              statisticsFile("") {}
        
        // Normalize move probabilities
        void normalizeProbs();
        
        // Validate configuration
        void validate() const;
    };
    
    // Constructor
    explicit GCMCModule(const Config& config = Config());
    ~GCMCModule();
    
    // Initialize with system state
    void initialize(model::montecarlo::MCState& state);
    
    // Fragment management
    int addFragmentType(const FragmentTemplate& tmpl);
    void setChemicalPotential(int typeId, double mu);
    void setActivity(int typeId, double activity);
    
    // Run simulation
    void runEquilibration(int steps = -1);
    void runProduction(int steps = -1);
    void runSteps(int nSteps);
    bool performMove();
    
    // Individual move types
    bool attemptInsertion(int typeId = -1);
    bool attemptDeletion(int typeId = -1);
    bool attemptTranslation();
    bool attemptRotation();
    bool attemptSwap();
    bool attemptRegrowth();
    bool attemptClusterMove();
    
    // Analysis
    const GCMCStats& getStatistics() const { return *statistics_; }
    void printStatistics() const;
    void saveStatistics(const std::string& filename) const;
    
    // Trajectory output
    void saveSnapshot();
    void saveTrajectory(const std::string& filename);
    
    // Energy calculation
    double calculateSystemEnergy();
    double calculateFragmentEnergy(int residueIdx);
    std::pair<double, double> calculateEnergyComponents();
    
    // Access to components
    GCMCEngine* getEngine() { return engine_.get(); }
    FragmentReservoir* getReservoir() { return reservoir_.get(); }
    CavityManager* getCavityManager() { return cavityManager_.get(); }
    
    // Configuration
    void setConfig(const Config& config) { config_ = config; }
    const Config& getConfig() const { return config_; }
    
    // Advanced features
    void enableAdaptiveBiasing();
    void setTargetDensity(double density);
    void enableFlatHistogram();
    
private:
    // Core components
    std::unique_ptr<GCMCEngine> engine_;
    std::unique_ptr<FragmentReservoir> reservoir_;
    std::unique_ptr<GCMCMoveSelector> moveSelector_;
    std::unique_ptr<GCMCBias> biasCalc_;
    std::unique_ptr<GCMCAcceptance> acceptCalc_;
    std::unique_ptr<GCMCStats> statistics_;
    
    // Bias components
    std::unique_ptr<CavityManager> cavityManager_;
    std::unique_ptr<ConfigBiasManager> configBias_;
    
    // State
    model::montecarlo::MCState* state_;
    Config config_;
    bool initialized_;
    int currentStep_;
    
    // Phase enum
    enum class Phase {
        EQUILIBRATION,
        PRODUCTION
    };
    Phase currentPhase_;
    
    // Helper methods
    void updateStatistics();
    void checkConvergence();
    void adaptBiasing();
};

/**
 * @brief Factory function to create configured GCMC module
 */
inline std::unique_ptr<GCMCModule> createGCMCModule(
    const GCMCModule::Config& config = GCMCModule::Config()) {
    return std::make_unique<GCMCModule>(config);
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_MODULE_HPP