#include "GCMCMoveSelector.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iomanip>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Constructor
GCMCMoveSelector::GCMCMoveSelector()
    : adaptiveProbabilities_(false),
      targetAcceptance_(0.3),
      adaptationRate_(0.01),
      adaptationInterval_(1000),
      movesSinceAdaptation_(0),
      rng_(0u),  // Deterministic initial value, will be set via setSeed()
      uniform_(0.0, 1.0) {

    // Initialize default probabilities
    probabilities_[MoveType::INSERT] = 0.25;
    probabilities_[MoveType::DELETE] = 0.25;
    probabilities_[MoveType::TRANSLATE] = 0.25;
    probabilities_[MoveType::ROTATE] = 0.25;
    probabilities_[MoveType::SWAP] = 0.0;
    probabilities_[MoveType::REGROWTH] = 0.0;
    probabilities_[MoveType::CLUSTER] = 0.0;
    probabilities_[MoveType::VOLUME_CHANGE] = 0.0;
    probabilities_[MoveType::IDENTITY_SWAP] = 0.0;

    updateCumulativeProbabilities();

    // Initialize statistics
    for (auto& [type, prob] : probabilities_) {
        statistics_[type] = MoveStats();
    }
}

// Destructor
GCMCMoveSelector::~GCMCMoveSelector() {
}

// Set move probabilities
void GCMCMoveSelector::setProbabilities(const std::map<MoveType, double>& probs) {
    probabilities_ = probs;
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Set single probability
void GCMCMoveSelector::setProbability(MoveType type, double prob) {
    probabilities_[type] = prob;
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Select next move type
GCMCMoveSelector::MoveType GCMCMoveSelector::selectMove() {
    double random = uniform_(rng_);
    
    for (const auto& [type, cumProb] : cumulativeProbabilities_) {
        if (random <= cumProb) {
            return type;
        }
    }
    
    // Fallback to translation (should never reach here)
    return MoveType::TRANSLATE;
}

// Select move weighted by acceptance rates
GCMCMoveSelector::MoveType GCMCMoveSelector::selectMoveWeighted() {
    // Calculate weights based on acceptance rates
    std::map<MoveType, double> weights;
    double totalWeight = 0.0;
    
    for (const auto& [type, prob] : probabilities_) {
        if (prob > 0) {
            double weight = calculateWeight(type);
            weights[type] = weight;
            totalWeight += weight;
        }
    }
    
    // Select based on weights
    double random = uniform_(rng_) * totalWeight;
    double cumWeight = 0.0;
    
    for (const auto& [type, weight] : weights) {
        cumWeight += weight;
        if (random <= cumWeight) {
            return type;
        }
    }
    
    return MoveType::TRANSLATE;
}

// Select fragment type for insertion/deletion
int GCMCMoveSelector::selectFragmentType(const std::vector<double>& chemicalPotentials) {
    if (chemicalPotentials.empty()) return -1;
    
    // Calculate probabilities based on chemical potentials
    std::vector<double> probs;
    double totalProb = 0.0;
    
    for (double mu : chemicalPotentials) {
        double prob = std::exp(mu / 8.314e-3 / 300.0);  // Simplified
        probs.push_back(prob);
        totalProb += prob;
    }
    
    // Select based on probabilities
    double random = uniform_(rng_) * totalProb;
    double cumProb = 0.0;
    
    for (size_t i = 0; i < probs.size(); ++i) {
        cumProb += probs[i];
        if (random <= cumProb) {
            return static_cast<int>(i);
        }
    }
    
    return 0;
}

// Select fragment type weighted by counts and activities
int GCMCMoveSelector::selectFragmentTypeWeighted(const std::vector<int>& counts,
                                                const std::vector<double>& activities) {
    if (counts.empty() || activities.empty()) return -1;
    
    // Calculate weights
    std::vector<double> weights;
    double totalWeight = 0.0;
    
    for (size_t i = 0; i < std::min(counts.size(), activities.size()); ++i) {
        double weight = activities[i] / (counts[i] + 1.0);
        weights.push_back(weight);
        totalWeight += weight;
    }
    
    // Select based on weights
    double random = uniform_(rng_) * totalWeight;
    double cumWeight = 0.0;
    
    for (size_t i = 0; i < weights.size(); ++i) {
        cumWeight += weights[i];
        if (random <= cumWeight) {
            return static_cast<int>(i);
        }
    }
    
    return 0;
}

// Select instance from available indices
int GCMCMoveSelector::selectInstance(const std::vector<int>& availableIndices) {
    if (availableIndices.empty()) return -1;
    
    std::uniform_int_distribution<int> dist(0, availableIndices.size() - 1);
    return availableIndices[dist(rng_)];
}

// Select swap pair
std::pair<int, int> GCMCMoveSelector::selectSwapPair(const std::vector<int>& type1Indices,
                                                     const std::vector<int>& type2Indices) {
    if (type1Indices.empty() || type2Indices.empty()) {
        return {-1, -1};
    }
    
    std::uniform_int_distribution<int> dist1(0, type1Indices.size() - 1);
    std::uniform_int_distribution<int> dist2(0, type2Indices.size() - 1);
    
    return {type1Indices[dist1(rng_)], type2Indices[dist2(rng_)]};
}

// Record attempt
void GCMCMoveSelector::recordAttempt(MoveType type) {
    statistics_[type].attempts++;
    movesSinceAdaptation_++;
    
    if (adaptiveProbabilities_ && movesSinceAdaptation_ >= adaptationInterval_) {
        updateAdaptiveProbabilities();
        movesSinceAdaptation_ = 0;
    }
}

// Record acceptance
void GCMCMoveSelector::recordAcceptance(MoveType type) {
    statistics_[type].accepted++;
}

// Record time
void GCMCMoveSelector::recordTime(MoveType type, double timeMs) {
    auto& stats = statistics_[type];
    stats.averageTime = (stats.averageTime * (stats.attempts - 1) + timeMs) / stats.attempts;
}

// Get statistics
const GCMCMoveSelector::MoveStats& GCMCMoveSelector::getStats(MoveType type) const {
    static const MoveStats emptyStats;
    auto it = statistics_.find(type);
    return (it != statistics_.end()) ? it->second : emptyStats;
}

// Print statistics
void GCMCMoveSelector::printStatistics() const {
    std::cout << "\n=== Move Statistics ===" << std::endl;
    std::cout << std::setw(15) << "Move Type" 
              << std::setw(10) << "Attempts"
              << std::setw(10) << "Accepted"
              << std::setw(12) << "Accept Rate"
              << std::setw(12) << "Avg Time(ms)" << std::endl;
    std::cout << std::string(59, '-') << std::endl;
    
    const char* moveNames[] = {
        "INSERT", "DELETE", "TRANSLATE", "ROTATE", "SWAP",
        "REGROWTH", "CLUSTER", "VOLUME_CHANGE", "IDENTITY_SWAP"
    };
    
    for (int i = 0; i < 9; ++i) {
        MoveType type = static_cast<MoveType>(i);
        auto it = statistics_.find(type);
        if (it != statistics_.end() && it->second.attempts > 0) {
            const auto& stats = it->second;
            std::cout << std::setw(15) << moveNames[i]
                      << std::setw(10) << stats.attempts
                      << std::setw(10) << stats.accepted
                      << std::setw(12) << std::fixed << std::setprecision(3) 
                      << stats.acceptanceRate()
                      << std::setw(12) << std::fixed << std::setprecision(2)
                      << stats.averageTime << std::endl;
        }
    }
}

// Reset statistics
void GCMCMoveSelector::resetStatistics() {
    for (auto& [type, stats] : statistics_) {
        stats = MoveStats();
    }
    movesSinceAdaptation_ = 0;
}

// Enable adaptive probabilities
void GCMCMoveSelector::enableAdaptiveProbabilities() {
    adaptiveProbabilities_ = true;
}

// Update adaptive probabilities
void GCMCMoveSelector::updateAdaptiveProbabilities() {
    // Adjust probabilities based on acceptance rates
    for (auto& [type, prob] : probabilities_) {
        if (prob > 0) {
            const auto& stats = statistics_[type];
            if (stats.attempts > 10) {  // Need enough statistics
                double acceptRate = stats.acceptanceRate();
                
                if (acceptRate < targetAcceptance_ - 0.1) {
                    // Too low - decrease probability
                    prob *= (1.0 - adaptationRate_);
                } else if (acceptRate > targetAcceptance_ + 0.1) {
                    // Too high - increase probability
                    prob *= (1.0 + adaptationRate_);
                }
            }
        }
    }
    
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Enable regrowth
void GCMCMoveSelector::enableRegrowth(double prob) {
    probabilities_[MoveType::REGROWTH] = prob;
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Enable cluster moves
void GCMCMoveSelector::enableClusterMoves(double prob) {
    probabilities_[MoveType::CLUSTER] = prob;
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Enable volume changes
void GCMCMoveSelector::enableVolumeChanges(double prob) {
    probabilities_[MoveType::VOLUME_CHANGE] = prob;
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Enable identity swaps
void GCMCMoveSelector::enableIdentitySwaps(double prob) {
    probabilities_[MoveType::IDENTITY_SWAP] = prob;
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Set insertion bias
void GCMCMoveSelector::setInsertionBias(int typeId, double bias) {
    insertionBiases_[typeId] = bias;
}

// Set deletion bias
void GCMCMoveSelector::setDeletionBias(int typeId, double bias) {
    deletionBiases_[typeId] = bias;
}

// Update cumulative probabilities
void GCMCMoveSelector::updateCumulativeProbabilities() {
    cumulativeProbabilities_.clear();
    double cumProb = 0.0;
    
    for (const auto& [type, prob] : probabilities_) {
        cumProb += prob;
        cumulativeProbabilities_[type] = cumProb;
    }
}

// Normalize probabilities
void GCMCMoveSelector::normalizeProbabilities() {
    double total = 0.0;
    for (const auto& [type, prob] : probabilities_) {
        total += prob;
    }
    
    if (total > 0) {
        for (auto& [type, prob] : probabilities_) {
            prob /= total;
        }
    }
}

// Calculate weight for move type
double GCMCMoveSelector::calculateWeight(MoveType type) const {
    const auto& stats = getStats(type);
    double baseProb = probabilities_.at(type);
    
    if (stats.attempts < 10) {
        // Not enough statistics - use base probability
        return baseProb;
    }
    
    // Weight by inverse of acceptance rate (encourage low-acceptance moves)
    double acceptRate = stats.acceptanceRate();
    double weight = baseProb;
    
    if (acceptRate < 0.1) {
        weight *= 2.0;  // Double weight for very low acceptance
    } else if (acceptRate > 0.5) {
        weight *= 0.5;  // Halve weight for very high acceptance
    }
    
    return weight;
}

// ============================================================================
// GCMCSmartMoveSelector Implementation
// ============================================================================

GCMCSmartMoveSelector::GCMCSmartMoveSelector()
    : GCMCMoveSelector(),
      learningRate_(0.1),
      discountFactor_(0.9),
      explorationRate_(0.1) {
}

// Learn from history
void GCMCSmartMoveSelector::learnFromHistory(const std::vector<MoveType>& history,
                                            const std::vector<bool>& accepted) {
    // Simple learning: increase probability of successful moves
    std::map<MoveType, int> successes;
    std::map<MoveType, int> attempts;
    
    for (size_t i = 0; i < history.size(); ++i) {
        attempts[history[i]]++;
        if (i < accepted.size() && accepted[i]) {
            successes[history[i]]++;
        }
    }
    
    // Update probabilities based on success rates
    for (auto& [type, count] : attempts) {
        if (count > 0) {
            double successRate = static_cast<double>(successes[type]) / count;
            double currentProb = probabilities_[type];
            probabilities_[type] = currentProb * (1 - learningRate_) + 
                                  successRate * learningRate_;
        }
    }
    
    normalizeProbabilities();
    updateCumulativeProbabilities();
}

// Predict best move
GCMCSmartMoveSelector::MoveType GCMCSmartMoveSelector::predictBestMove(
    double currentEnergy, int nMolecules, double temperature) {
    
    // Suppress unused parameter warning
    (void)temperature;
    
    int state = discretizeState(currentEnergy, nMolecules);
    
    // Epsilon-greedy selection
    if (uniform_(rng_) < explorationRate_) {
        // Explore: random move
        return selectMove();
    }
    
    // Exploit: select best Q-value
    MoveType bestMove = MoveType::TRANSLATE;
    double bestQ = -1e10;
    
    for (int i = 0; i < 9; ++i) {
        MoveType type = static_cast<MoveType>(i);
        auto key = std::make_pair(state, type);
        
        if (qValues_.find(key) != qValues_.end()) {
            if (qValues_[key] > bestQ) {
                bestQ = qValues_[key];
                bestMove = type;
            }
        }
    }
    
    return bestMove;
}

// Update Q-values
void GCMCSmartMoveSelector::updateQValues(MoveType type, double reward) {
    // Simplified Q-learning update
    // In practice, would need previous state and action
    int state = 0;  // Placeholder
    auto key = std::make_pair(state, type);
    
    double oldQ = (qValues_.find(key) != qValues_.end()) ? qValues_[key] : 0.0;
    double maxNextQ = 0.0;  // Would calculate max Q over next actions
    
    qValues_[key] = oldQ + learningRate_ * (reward + discountFactor_ * maxNextQ - oldQ);
}

// Discretize state for Q-learning
int GCMCSmartMoveSelector::discretizeState(double energy, int nMolecules) {
    // Simple discretization
    int energyBin = static_cast<int>(energy / 100.0);
    int moleculeBin = nMolecules / 10;
    
    return energyBin * 100 + moleculeBin;
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc