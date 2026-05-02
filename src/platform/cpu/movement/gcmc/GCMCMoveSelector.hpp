#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_MOVE_SELECTOR_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_MOVE_SELECTOR_HPP

#include <random>
#include <vector>
#include <map>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Selects move types based on configured probabilities
 *
 * This class manages the selection of different move types in GCMC
 * simulations based on user-defined probabilities. It supports:
 * - Standard moves (insert, delete, translate, rotate)
 * - Advanced moves (swap, regrowth, cluster)
 * - Adaptive probability adjustment
 */
class GCMCMoveSelector {
public:
    // Move types
    enum class MoveType {
        INSERT,
        DELETE,
        TRANSLATE,
        ROTATE,
        SWAP,
        REGROWTH,
        CLUSTER,
        VOLUME_CHANGE,  // For NPT-GCMC
        IDENTITY_SWAP   // For mixture simulations
    };

    // Move statistics
    struct MoveStats {
        int attempts = 0;
        int accepted = 0;
        double averageTime = 0.0;  // ms

        double acceptanceRate() const {
            return attempts > 0 ?
                   static_cast<double>(accepted) / attempts : 0.0;
        }
    };

    // Constructor
    GCMCMoveSelector();
    ~GCMCMoveSelector();

    // Set move probabilities (must sum to 1.0)
    void setProbabilities(const std::map<MoveType, double>& probs);
    void setProbability(MoveType type, double prob);

    // Select next move type
    MoveType selectMove();
    MoveType selectMoveWeighted();  // Weight by acceptance rates

    // Fragment type selection for insertion/deletion
    int selectFragmentType(const std::vector<double>& chemicalPotentials);
    int selectFragmentTypeWeighted(const std::vector<int>& counts,
                                   const std::vector<double>& activities);

    // Instance selection for moves
    int selectInstance(const std::vector<int>& availableIndices);
    std::pair<int, int> selectSwapPair(const std::vector<int>& type1Indices,
                                       const std::vector<int>& type2Indices);

    // Update statistics
    void recordAttempt(MoveType type);
    void recordAcceptance(MoveType type);
    void recordTime(MoveType type, double timeMs);

    // Get statistics
    const MoveStats& getStats(MoveType type) const;
    void printStatistics() const;
    void resetStatistics();

    // Adaptive features
    void enableAdaptiveProbabilities();
    void updateAdaptiveProbabilities();
    void setTargetAcceptance(double target) { targetAcceptance_ = target; }

    // Special move configurations
    void enableRegrowth(double prob);
    void enableClusterMoves(double prob);
    void enableVolumeChanges(double prob);
    void enableIdentitySwaps(double prob);

    // Bias for specific fragment types
    void setInsertionBias(int typeId, double bias);
    void setDeletionBias(int typeId, double bias);

    // Random seed setting
    void setSeed(unsigned int seed) { rng_.seed(seed); }

protected:
    // Probabilities
    std::map<MoveType, double> probabilities_;
    std::map<MoveType, double> cumulativeProbabilities_;

    // Statistics
    std::map<MoveType, MoveStats> statistics_;

    // Fragment-specific biases
    std::map<int, double> insertionBiases_;
    std::map<int, double> deletionBiases_;

    // Adaptive parameters
    bool adaptiveProbabilities_;
    double targetAcceptance_;
    double adaptationRate_;
    int adaptationInterval_;
    int movesSinceAdaptation_;

    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;

    // Helper methods
    void updateCumulativeProbabilities();
    void normalizeProbabilities();
    double calculateWeight(MoveType type) const;
};

/**
 * @brief Advanced move selector with machine learning capabilities
 */
class GCMCSmartMoveSelector : public GCMCMoveSelector {
public:
    GCMCSmartMoveSelector();

    // Learn from history
    void learnFromHistory(const std::vector<MoveType>& history,
                          const std::vector<bool>& accepted);

    // Predict best move type
    MoveType predictBestMove(double currentEnergy,
                             int nMolecules,
                             double temperature);

    // Reinforcement learning
    void updateQValues(MoveType type, double reward);

private:
    // Q-learning parameters
    std::map<std::pair<int, MoveType>, double> qValues_;  // (state, action) -> Q
    double learningRate_;
    double discountFactor_;
    double explorationRate_;

    // State discretization
    int discretizeState(double energy, int nMolecules);
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_MOVE_SELECTOR_HPP
