#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_CONFIG_BIAS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_CONFIG_BIAS_HPP

#include <vector>
#include <memory>
#include <algorithm>
#include "../common/MovementUtils.hpp"

// Forward declarations
namespace pygcmc {
namespace model {
namespace montecarlo {
    class MCState;
    class MCResidue;
}
}
}

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCState;
using model::montecarlo::MCResidue;

// Energy interface placeholder
class EnergyInterface {
public:
    virtual ~EnergyInterface() = default;
    virtual double calculateEnergy(const MCState& state) = 0;
};

/**
 * Configuration for a trial orientation/position
 */
struct Configuration {
    Quaternion rotation;         // Rotation quaternion
    Vector3 translation;         // Translation vector
    double energy;              // Energy of this configuration
    double probability;         // Probability weight
    int index;                  // Configuration index
    
    Configuration() : energy(0.0), probability(0.0), index(0) {}
};

/**
 * Configurational Bias Manager for improved acceptance rates
 * Generates multiple trial configurations and selects based on Boltzmann weights
 */
class ConfigBiasManager {
public:
    // Constructor
    explicit ConfigBiasManager(int numTrials = 10);
    
    // Destructor
    ~ConfigBiasManager();
    
    // Main configuration selection function
    Configuration selectConfiguration(
        MCState& state,
        int residueIdx,
        int numTrials,
        EnergyInterface* energyCalc,
        double beta
    );
    
    // Generate trial configurations for rotation
    std::vector<Configuration> generateRotationConfigurations(
        const MCState& state,
        int residueIdx,
        int numTrials
    );
    
    // Generate trial configurations for insertion
    std::vector<Configuration> generateInsertionConfigurations(
        const Vector3& basePosition,
        int numTrials,
        double translationRange
    );
    
    // Calculate Boltzmann probabilities for configurations
    void calculateProbabilities(
        std::vector<Configuration>& configs,
        double beta,
        bool useLogSpace = true
    );
    
    // Select configuration based on probability distribution
    int selectByProbability(const std::vector<Configuration>& configs);
    
    // Calculate bias correction factor
    double calculateBiasFactor(
        const std::vector<Configuration>& configs,
        int selectedIndex
    );
    
    // Configuration
    void setNumTrials(int trials) { numTrials_ = trials; }
    int getNumTrials() const { return numTrials_; }
    
    void setTranslationRange(double range) { translationRange_ = range; }
    double getTranslationRange() const { return translationRange_; }
    
    void setIncludeTranslation(bool include) { includeTranslation_ = include; }
    bool getIncludeTranslation() const { return includeTranslation_; }
    
    // Statistics
    struct Statistics {
        int totalSelections = 0;
        double averageNumConfigs = 0.0;
        double averageBiasFactor = 0.0;
        double averageAcceptance = 0.0;
        int minEnergySelected = 0;  // Times minimum energy config was selected
        std::vector<double> energyDistribution;
        std::vector<double> probabilityDistribution;
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
    // Advanced features
    struct ConfigurationSet {
        std::vector<Configuration> configs;
        int bestIndex = -1;         // Index of lowest energy
        double minEnergy = 0.0;      // Minimum energy
        double partitionFunction = 0.0;  // Sum of Boltzmann weights
    };
    
    // Generate and evaluate a full set of configurations
    ConfigurationSet generateConfigurationSet(
        MCState& state,
        int residueIdx,
        EnergyInterface* energyCalc,
        double beta
    );
    
    // Apply configuration to state (temporarily)
    void applyConfiguration(
        MCState& state,
        int residueIdx,
        const Configuration& config
    );
    
    // Restore original configuration
    void restoreConfiguration(
        MCState& state,
        int residueIdx,
        const Configuration& original
    );
    
private:
    // Configuration parameters
    int numTrials_;                  // Number of trial configurations
    double translationRange_;        // Range for translation in Angstroms
    bool includeTranslation_;       // Include translation in config bias
    
    // Statistics
    mutable Statistics stats_;
    
    // Helper functions
    void normalizeConfigurations(std::vector<Configuration>& configs);
    double computePartitionFunction(const std::vector<Configuration>& configs, double beta);
    void updateStatistics(const ConfigurationSet& configSet, int selectedIndex);
    
    // Random configuration generation
    Quaternion generateRandomRotation();
    Vector3 generateRandomTranslation(double range);
    
    // Energy evaluation helpers
    double evaluateConfiguration(
        MCState& state,
        int residueIdx,
        const Configuration& config,
        EnergyInterface* energyCalc
    );
};

/**
 * Specialized configurational bias for rotation moves
 */
class ConfigBiasRotation {
public:
    // Constructor
    explicit ConfigBiasRotation(ConfigBiasManager* manager);
    
    // Perform rotation with configurational bias
    bool performRotation(
        MCState& state,
        int residueIdx,
        EnergyInterface* energyCalc,
        const MovementParams& params
    );
    
    // Calculate acceptance probability with bias correction
    double calculateAcceptanceProbability(
        double deltaE,
        double beta,
        double biasFactor
    );
    
    // Statistics
    struct Statistics {
        int attempts = 0;
        int accepts = 0;
        double averageBiasFactor = 0.0;
        double acceptanceRate() const {
            return attempts > 0 ? static_cast<double>(accepts) / attempts : 0.0;
        }
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
private:
    ConfigBiasManager* manager_;
    Statistics stats_;
    
    // Helper functions
    void saveOriginalConfiguration(MCState& state, int residueIdx, Configuration& config);
    void updateStatistics(bool accepted, double biasFactor);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_CONFIG_BIAS_HPP