#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ACCEPTANCE_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ACCEPTANCE_HPP

#include <cmath>
#include <random>
#include <map>
#include <vector>
#include "../common/GrandCanonicalTerms.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Calculates acceptance probabilities for GCMC moves
 * 
 * This class implements various acceptance criteria:
 * - Metropolis criterion
 * - Grand canonical acceptance
 * - NPT-GCMC acceptance
 * - Multi-component acceptance
 */
class GCMCAcceptance {
public:
    enum class MoveType {
        INSERTION,
        DELETION
    };

    struct GrandCanonicalInsertionTerms {
        int typeId = -1;
        int countBefore = 0;
        double deltaE = 0.0;
        double cavityFraction = 1.0;
        double lambdaNm = 1.0;
        double rosenbluthWeight = 1.0;
        int cbmcTrials = 1;
        double proposalLogRatio = 0.0; // log(p_fwd) - log(p_rev)
    };

    struct GrandCanonicalDeletionTerms {
        int typeId = -1;
        int countBefore = 0;
        double deltaE = 0.0;
        double cavityFraction = 1.0;
        double lambdaNm = 1.0;
        double rosenbluthWeight = 1.0;
        int cbmcTrials = 1;
        double proposalLogRatio = 0.0; // log(p_fwd) - log(p_rev)
    };

    // Acceptance criteria types
    enum class CriterionType {
        METROPOLIS,           // Standard Metropolis
        GRAND_CANONICAL,      // μVT ensemble
        ISOTHERMAL_ISOBARIC, // NPT ensemble
        GIBBS_ENSEMBLE,      // Gibbs ensemble
        WANG_LANDAU,         // Flat histogram
        TRANSITION_MATRIX    // Transition matrix MC
    };
    
    // Constructor
    GCMCAcceptance();
    ~GCMCAcceptance();
    
    // Set ensemble parameters
    void setTemperature(double T) { temperature_ = T; }
    void setPressure(double P) { pressure_ = P; }
    void setVolume(double V) { volume_ = V; }
    double getVolume() const { return volume_; }  // Get current volume
    void setChemicalPotential(int typeId, double mu);
    void setActivity(int typeId, double activity);
    void setThermalLambda(int typeId, double lambdaNm);
    double getThermalLambda(int typeId) const;
    double getActivity(int typeId) const;
    
    // Calculate acceptance probability
    double calculateInsertionProbability(
        int typeId,
        int currentNumber,
        double deltaE,
        double bias = 1.0
    );
    
    double calculateDeletionProbability(
        int typeId,
        int currentNumber,
        double deltaE,
        double bias = 1.0
    );

    double calculateInsertionProbabilityDetailed(
        const GrandCanonicalInsertionTerms& terms,
        double* logRatioOut = nullptr
    );

    double calculateDeletionProbabilityDetailed(
        const GrandCanonicalDeletionTerms& terms,
        double* logRatioOut = nullptr
    );

    // Unified evaluator shared by MovementModule and CLI engine
    static gcmc::GrandCanonicalEvaluation evaluate(
        const gcmc::GrandCanonicalTerms& terms,
        MoveType moveType
    );
    
    double calculateTranslationProbability(
        double deltaE,
        double bias = 1.0
    );
    
    double calculateRotationProbability(
        double deltaE,
        double bias = 1.0
    );
    
    double calculateSwapProbability(
        int type1, int type2,
        int n1, int n2,
        double deltaE,
        double bias = 1.0
    );
    
    double calculateVolumeChangeProbability(
        double oldVolume,
        double newVolume,
        int nMolecules,
        double deltaE
    );
    
    // General acceptance calculation
    double calculateAcceptance(
        CriterionType criterion,
        double deltaE,
        double bias = 1.0,
        double additionalFactor = 1.0
    );
    
    // Accept or reject based on probability
    bool acceptMove(double probability);
    bool acceptMetropolis(double deltaE, double bias = 1.0);
    
    // Advanced acceptance criteria
    double calculateGibbsAcceptance(
        int boxFrom, int boxTo,
        int typeId,
        double deltaE
    );
    
    double calculateWangLandauAcceptance(
        double currentBias,
        double newBias
    );
    
    double calculateTransitionMatrixAcceptance(
        int oldState,
        int newState,
        const std::vector<std::vector<double>>& transitionMatrix
    );
    
    // Detailed balance checks
    bool checkDetailedBalance(
        double forwardProb,
        double reverseProb,
        double tolerance = 1e-6
    );
    
    // Configuration
    void setCriterion(CriterionType type) { criterionType_ = type; }
    void setSeed(unsigned int seed) { rng_.seed(seed); }
    
    // Statistics
    int getTotalDecisions() const { return totalDecisions_; }
    int getAcceptedMoves() const { return acceptedMoves_; }
    double getAverageAcceptance() const {
        return totalDecisions_ > 0 ? 
               static_cast<double>(acceptedMoves_) / totalDecisions_ : 0.0;
    }
    
protected:
    // Ensemble parameters
    double temperature_;
    double pressure_;
    double volume_;
    std::map<int, double> chemicalPotentials_;
    std::map<int, double> activities_;
    std::map<int, double> thermalLambdaNm_;
    
    // Configuration
    CriterionType criterionType_;
    
    // Constants
    static constexpr double kB = 8.314e-3;  // kJ/(mol·K)
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;
    
    // Statistics
    int totalDecisions_;
    int acceptedMoves_;
    
    // Helper methods
    double getBeta() const { return 1.0 / (kB * temperature_); }
    
private:
    double getIdealGasContribution(int n, double V) const;
    double getDeBroglieWavelength(double mass) const;
    double safeLog(double value) const;
    double safeExp(double logValue) const;
    double getChemicalPotentialInternal(int typeId) const;
};

/**
 * @brief Smart acceptance calculator with machine learning
 */
class GCMCSmartAcceptance : public GCMCAcceptance {
public:
    GCMCSmartAcceptance();
    
    // Learn optimal acceptance from history
    void learnFromHistory(
        const std::vector<double>& deltaEs,
        const std::vector<bool>& accepted,
        const std::vector<double>& systemProperties
    );
    
    // Predict acceptance probability
    double predictAcceptance(
        double deltaE,
        const std::vector<double>& systemProperties
    );
    
    // Adjust criterion dynamically
    void adjustCriterion(double targetAcceptance);
    
private:
    // Machine learning model parameters
    std::vector<double> weights_;
    double learningRate_;
    
    // Feature extraction
    std::vector<double> extractFeatures(
        double deltaE,
        const std::vector<double>& properties
    );
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ACCEPTANCE_HPP
