#include "ConfigBias.hpp"
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../../../simulation/simulation.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

ConfigBiasManager::ConfigBiasManager(int numTrials)
    : numTrials_(numTrials),
      translationRange_(0.5),
      includeTranslation_(false) {
    resetStatistics();
}

ConfigBiasManager::~ConfigBiasManager() = default;

Configuration ConfigBiasManager::selectConfiguration(
    MCState& state,
    int residueIdx,
    int numTrials,
    EnergyInterface* energyCalc,
    double beta) {
    
    // Generate configurations
    auto configs = generateRotationConfigurations(state, residueIdx, numTrials);
    
    // Evaluate energies (would need energy calculator)
    for (auto& config : configs) {
        config.energy = evaluateConfiguration(state, residueIdx, config, energyCalc);
    }
    
    // Calculate probabilities
    calculateProbabilities(configs, beta, true);
    
    // Select based on probability
    int selectedIdx = selectByProbability(configs);
    
    if (selectedIdx >= 0 && selectedIdx < static_cast<int>(configs.size())) {
        return configs[selectedIdx];
    }
    
    // Return default if selection failed
    return Configuration();
}

std::vector<Configuration> ConfigBiasManager::generateRotationConfigurations(
    const MCState& /*state*/,
    int /*residueIdx*/,
    int numTrials) {
    
    std::vector<Configuration> configs;
    configs.reserve(numTrials);
    
    for (int i = 0; i < numTrials; ++i) {
        Configuration config;
        config.index = i;
        config.rotation = generateRandomRotation();
        
        if (includeTranslation_) {
            config.translation = generateRandomTranslation(translationRange_);
        }
        
        configs.push_back(config);
    }
    
    return configs;
}

std::vector<Configuration> ConfigBiasManager::generateInsertionConfigurations(
    const Vector3& basePosition,
    int numTrials,
    double translationRange) {
    
    std::vector<Configuration> configs;
    configs.reserve(numTrials);
    
    for (int i = 0; i < numTrials; ++i) {
        Configuration config;
        config.index = i;
        config.rotation = generateRandomRotation();
        config.translation = basePosition + generateRandomTranslation(translationRange);
        configs.push_back(config);
    }
    
    return configs;
}

void ConfigBiasManager::calculateProbabilities(
    std::vector<Configuration>& configs,
    double beta,
    bool useLogSpace) {
    
    if (configs.empty()) return;
    
    if (useLogSpace) {
        // Find minimum energy for numerical stability
        double minEnergy = std::numeric_limits<double>::max();
        for (const auto& config : configs) {
            minEnergy = std::min(minEnergy, config.energy);
        }
        
        // Calculate log probabilities
        std::vector<double> logProbs;
        logProbs.reserve(configs.size());
        
        for (const auto& config : configs) {
            logProbs.push_back(-beta * (config.energy - minEnergy));
        }
        
        // Calculate log sum
        double logSum = utils::LogSpaceCalculator::logSumExp(logProbs);
        
        // Normalize probabilities
        for (size_t i = 0; i < configs.size(); ++i) {
            configs[i].probability = std::exp(logProbs[i] - logSum);
        }
    } else {
        // Direct calculation
        double sum = 0.0;
        
        // Find minimum energy for better numerical stability
        double minEnergy = std::numeric_limits<double>::max();
        for (const auto& config : configs) {
            minEnergy = std::min(minEnergy, config.energy);
        }
        
        // Calculate unnormalized probabilities
        for (auto& config : configs) {
            config.probability = std::exp(-beta * (config.energy - minEnergy));
            sum += config.probability;
        }
        
        // Normalize
        if (sum > 0.0) {
            for (auto& config : configs) {
                config.probability /= sum;
            }
        }
    }
    
    // Update statistics
    stats_.totalSelections++;
    stats_.averageNumConfigs = (stats_.averageNumConfigs * (stats_.totalSelections - 1) + configs.size()) / 
                               stats_.totalSelections;
}

int ConfigBiasManager::selectByProbability(const std::vector<Configuration>& configs) {
    if (configs.empty()) return -1;
    
    // Generate random number
    double r = utils::RandomUtils::uniform(0.0, 1.0);
    
    // Select based on cumulative probability
    double cumsum = 0.0;
    for (size_t i = 0; i < configs.size(); ++i) {
        cumsum += configs[i].probability;
        if (r <= cumsum) {
            // Check if this is the minimum energy configuration
            double minEnergy = std::numeric_limits<double>::max();
            for (const auto& config : configs) {
                minEnergy = std::min(minEnergy, config.energy);
            }
            if (std::abs(configs[i].energy - minEnergy) < 1e-10) {
                stats_.minEnergySelected++;
            }
            return static_cast<int>(i);
        }
    }
    
    // Return last if cumsum rounding issue
    return static_cast<int>(configs.size()) - 1;
}

double ConfigBiasManager::calculateBiasFactor(
    const std::vector<Configuration>& configs,
    int selectedIndex) {
    
    if (selectedIndex < 0 || selectedIndex >= static_cast<int>(configs.size())) {
        return 1.0;
    }
    
    // Rosenbluth factor = n_trials * P(selected)
    double biasFactor = configs.size() * configs[selectedIndex].probability;
    
    // Update statistics
    stats_.averageBiasFactor = (stats_.averageBiasFactor * (stats_.totalSelections - 1) + biasFactor) / 
                              stats_.totalSelections;
    
    // Ensure non-zero
    if (biasFactor < 1e-10) {
        biasFactor = 1e-10;
    }
    
    return biasFactor;
}

ConfigBiasManager::ConfigurationSet ConfigBiasManager::generateConfigurationSet(
    MCState& state,
    int residueIdx,
    EnergyInterface* energyCalc,
    double beta) {
    
    ConfigurationSet configSet;
    
    // Generate configurations
    configSet.configs = generateRotationConfigurations(state, residueIdx, numTrials_);
    
    // Evaluate energies
    configSet.minEnergy = std::numeric_limits<double>::max();
    for (auto& config : configSet.configs) {
        config.energy = evaluateConfiguration(state, residueIdx, config, energyCalc);
        if (config.energy < configSet.minEnergy) {
            configSet.minEnergy = config.energy;
            configSet.bestIndex = config.index;
        }
    }
    
    // Calculate probabilities
    calculateProbabilities(configSet.configs, beta, true);
    
    // Calculate partition function
    configSet.partitionFunction = computePartitionFunction(configSet.configs, beta);
    
    return configSet;
}

void ConfigBiasManager::resetStatistics() {
    stats_ = Statistics();
    stats_.energyDistribution.clear();
    stats_.probabilityDistribution.clear();
}

// Private helper functions

void ConfigBiasManager::normalizeConfigurations(std::vector<Configuration>& configs) {
    double sum = 0.0;
    for (const auto& config : configs) {
        sum += config.probability;
    }
    
    if (sum > 0.0) {
        for (auto& config : configs) {
            config.probability /= sum;
        }
    }
}

double ConfigBiasManager::computePartitionFunction(const std::vector<Configuration>& configs, double beta) {
    double Z = 0.0;
    
    // Find minimum energy for numerical stability
    double minEnergy = std::numeric_limits<double>::max();
    for (const auto& config : configs) {
        minEnergy = std::min(minEnergy, config.energy);
    }
    
    // Sum Boltzmann weights
    for (const auto& config : configs) {
        Z += std::exp(-beta * (config.energy - minEnergy));
    }
    
    // Multiply back the factor
    Z *= std::exp(-beta * minEnergy);
    
    return Z;
}

void ConfigBiasManager::updateStatistics(const ConfigurationSet& configSet, int selectedIndex) {
    // Record energy distribution
    for (const auto& config : configSet.configs) {
        stats_.energyDistribution.push_back(config.energy);
        stats_.probabilityDistribution.push_back(config.probability);
    }
    
    // Update averages
    if (selectedIndex >= 0 && selectedIndex < static_cast<int>(configSet.configs.size())) {
        double acceptance = configSet.configs[selectedIndex].probability;
        stats_.averageAcceptance = (stats_.averageAcceptance * (stats_.totalSelections - 1) + acceptance) / 
                                  stats_.totalSelections;
    }
}

Quaternion ConfigBiasManager::generateRandomRotation() {
    return utils::RotationUtils::generateRandomQuaternion();
}

Vector3 ConfigBiasManager::generateRandomTranslation(double range) {
    return utils::RandomUtils::randomVector(range);
}

double ConfigBiasManager::evaluateConfiguration(
    MCState& state,
    int /*residueIdx*/,
    const Configuration& /*config*/,
    EnergyInterface* /*energyCalc*/) {
    
    // This is a simplified version - actual implementation would:
    // 1. Save current configuration
    // 2. Apply the configuration (rotation/translation)
    // 3. Calculate energy
    // 4. Restore original configuration
    
    // For now, use direct energy calculation
    simulation::Simulation::computeSystemEnergyCutoff(state);
    
    double energy = 0.0;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        energy += state.residues[i].energy_vdw;
        energy += state.residues[i].energy_elec;
    }
    
    return energy * 0.5;  // Account for double counting
}

void ConfigBiasManager::applyConfiguration(
    MCState& /*state*/,
    int /*residueIdx*/,
    const Configuration& /*config*/) {
    
    // Apply rotation and translation to residue
    // Implementation would manipulate atoms in state.atoms[]
    // This is handled in the specific move classes
}

void ConfigBiasManager::restoreConfiguration(
    MCState& /*state*/,
    int /*residueIdx*/,
    const Configuration& /*original*/) {
    
    // Restore original configuration
    // Implementation would restore atoms in state.atoms[]
}

// ConfigBiasRotation implementation

ConfigBiasRotation::ConfigBiasRotation(ConfigBiasManager* manager)
    : manager_(manager) {
    resetStatistics();
}

bool ConfigBiasRotation::performRotation(
    MCState& state,
    int residueIdx,
    EnergyInterface* energyCalc,
    const MovementParams& params) {
    
    // Generate configuration set
    auto configSet = manager_->generateConfigurationSet(state, residueIdx, energyCalc, params.beta);
    
    // Select configuration
    int selectedIdx = manager_->selectByProbability(configSet.configs);
    
    if (selectedIdx < 0) {
        return false;
    }
    
    // Calculate bias factor
    double biasFactor = manager_->calculateBiasFactor(configSet.configs, selectedIdx);
    
    // Apply selected configuration
    manager_->applyConfiguration(state, residueIdx, configSet.configs[selectedIdx]);
    
    // Calculate energy change
    double deltaE = configSet.configs[selectedIdx].energy - configSet.minEnergy;
    
    // Calculate acceptance probability with bias correction
    double acceptProb = calculateAcceptanceProbability(deltaE, params.beta, biasFactor);
    
    // Accept or reject
    bool accepted = utils::RandomUtils::metropolisAccept(acceptProb);
    
    if (!accepted) {
        // Restore original configuration
        Configuration original;
        saveOriginalConfiguration(state, residueIdx, original);
        manager_->restoreConfiguration(state, residueIdx, original);
    }
    
    // Update statistics
    updateStatistics(accepted, biasFactor);
    
    return accepted;
}

double ConfigBiasRotation::calculateAcceptanceProbability(
    double deltaE,
    double beta,
    double biasFactor) {
    
    // Bias-corrected acceptance probability
    // A = min(1, exp(-β*ΔE) / bias_factor)
    return std::min(1.0, std::exp(-beta * deltaE) / biasFactor);
}

void ConfigBiasRotation::resetStatistics() {
    stats_ = Statistics();
}

void ConfigBiasRotation::saveOriginalConfiguration(MCState& /*state*/, int /*residueIdx*/, Configuration& /*config*/) {
    // Save current configuration
    // This would save atom positions for the residue
    // Implementation handled in specific move classes
}

void ConfigBiasRotation::updateStatistics(bool accepted, double biasFactor) {
    stats_.attempts++;
    if (accepted) {
        stats_.accepts++;
    }
    
    stats_.averageBiasFactor = (stats_.averageBiasFactor * (stats_.attempts - 1) + biasFactor) / stats_.attempts;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc