#include "GCMCAcceptance.hpp"
#include <algorithm>
#include <numeric>
#include <limits>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Constructor
GCMCAcceptance::GCMCAcceptance()
    : temperature_(300.0),
      pressure_(1.0),
      volume_(1000.0),
      criterionType_(CriterionType::GRAND_CANONICAL),
      rng_(0u),  // Deterministic initial value, will be set via setSeed()
      uniform_(0.0, 1.0),
      totalDecisions_(0),
      acceptedMoves_(0) {
}

// Destructor
GCMCAcceptance::~GCMCAcceptance() {
}

// Set chemical potential for a type
void GCMCAcceptance::setChemicalPotential(int typeId, double mu) {
    chemicalPotentials_[typeId] = mu;
    // Update activity
    double beta = getBeta();
    activities_[typeId] = std::exp(beta * mu);
    if (thermalLambdaNm_.find(typeId) == thermalLambdaNm_.end()) {
        thermalLambdaNm_[typeId] = 1.0;
    }
}

// Set activity for a type
void GCMCAcceptance::setActivity(int typeId, double activity) {
    activities_[typeId] = activity;
    // Update chemical potential
    double beta = getBeta();
    chemicalPotentials_[typeId] = std::log(activity) / beta;
    if (thermalLambdaNm_.find(typeId) == thermalLambdaNm_.end()) {
        thermalLambdaNm_[typeId] = 1.0;
    }
}

void GCMCAcceptance::setThermalLambda(int typeId, double lambdaNm) {
    if (lambdaNm > 0.0) {
        thermalLambdaNm_[typeId] = lambdaNm;
    } else {
        thermalLambdaNm_[typeId] = 1.0;
    }
}

double GCMCAcceptance::getThermalLambda(int typeId) const {
    auto it = thermalLambdaNm_.find(typeId);
    if (it != thermalLambdaNm_.end() && it->second > 0.0) {
        return it->second;
    }
    return 1.0;
}

double GCMCAcceptance::getActivity(int typeId) const {
    auto it = activities_.find(typeId);
    if (it != activities_.end()) {
        return it->second;
    }
    auto muIt = chemicalPotentials_.find(typeId);
    if (muIt != chemicalPotentials_.end()) {
        double beta = getBeta();
        return std::exp(beta * muIt->second);
    }
    return 1.0;
}

// Calculate insertion probability
double GCMCAcceptance::calculateInsertionProbability(
    int typeId,
    int currentNumber,
    double deltaE,
    double bias) {
    
    double beta = getBeta();
    double activity = activities_[typeId];
    
    // Grand canonical acceptance: min(1, (zV/(N+1)) * exp(-beta*deltaE) * bias)
    double prefactor = activity * volume_ / (currentNumber + 1);
    double boltzmann = std::exp(-beta * deltaE);
    
    return std::min(1.0, prefactor * boltzmann * bias);
}

double GCMCAcceptance::calculateInsertionProbabilityDetailed(
    const GrandCanonicalInsertionTerms& terms,
    double* logRatioOut) {
    
    if (terms.typeId < 0) {
        if (logRatioOut) {
            *logRatioOut = -std::numeric_limits<double>::infinity();
        }
        return 0.0;
    }

    gcmc::GrandCanonicalTerms unified;
    unified.speciesId = terms.typeId;
    unified.countBefore = terms.countBefore;
    unified.countAfter = terms.countBefore + 1;
    unified.beta = getBeta();
    unified.chemicalPotential = getChemicalPotentialInternal(terms.typeId);
    unified.deltaEnergy = terms.deltaE;

    double volume = volume_ > 0.0 ? volume_ : 1.0;
    unified.logVolume = safeLog(volume);

    double cavityFraction = std::max(terms.cavityFraction, 1e-30);
    unified.logCavityForward = safeLog(cavityFraction);
    unified.logCavityReverse = 0.0;

    double lambdaNm = terms.lambdaNm > 0.0 ? terms.lambdaNm : getThermalLambda(terms.typeId);
    unified.logLambda3 = 3.0 * safeLog(lambdaNm);

    unified.logProposalForward = 0.0;
    unified.logProposalReverse = terms.proposalLogRatio;

    double rosen = std::max(terms.rosenbluthWeight, 1e-30);
    unified.logRosenbluthForward = safeLog(rosen);
    unified.logRosenbluthReverse = 0.0;

    unified.logExtraForward = 0.0;
    unified.logExtraReverse = 0.0;

    auto eval = evaluate(unified, MoveType::INSERTION);
    if (logRatioOut) {
        *logRatioOut = eval.logRatio;
    }
    return eval.probability;
}

// Calculate deletion probability
double GCMCAcceptance::calculateDeletionProbability(
    int typeId,
    int currentNumber,
    double deltaE,
    double bias) {
    
    if (currentNumber == 0) return 0.0;
    
    double beta = getBeta();
    double activity = activities_[typeId];
    
    // Grand canonical acceptance: min(1, (N/(zV)) * exp(-beta*deltaE) * bias)
    double prefactor = currentNumber / (activity * volume_);
    double boltzmann = std::exp(-beta * deltaE);
    
    return std::min(1.0, prefactor * boltzmann * bias);
}

double GCMCAcceptance::calculateDeletionProbabilityDetailed(
    const GrandCanonicalDeletionTerms& terms,
    double* logRatioOut) {
    
    if (terms.typeId < 0 || terms.countBefore <= 0) {
        if (logRatioOut) {
            *logRatioOut = -std::numeric_limits<double>::infinity();
        }
        return 0.0;
    }

    gcmc::GrandCanonicalTerms unified;
    unified.speciesId = terms.typeId;
    unified.countBefore = terms.countBefore;
    unified.countAfter = std::max(terms.countBefore - 1, 0);
    unified.beta = getBeta();
    unified.chemicalPotential = getChemicalPotentialInternal(terms.typeId);
    unified.deltaEnergy = terms.deltaE;

    double volume = volume_ > 0.0 ? volume_ : 1.0;
    unified.logVolume = safeLog(volume);

    double cavityFraction = std::max(terms.cavityFraction, 1e-30);
    unified.logCavityForward = 0.0;
    unified.logCavityReverse = safeLog(cavityFraction);

    double lambdaNm = terms.lambdaNm > 0.0 ? terms.lambdaNm : getThermalLambda(terms.typeId);
    unified.logLambda3 = 3.0 * safeLog(lambdaNm);

    unified.logProposalForward = 0.0;
    unified.logProposalReverse = terms.proposalLogRatio;

    double rosen = std::max(terms.rosenbluthWeight, 1e-30);
    unified.logRosenbluthForward = 0.0;
    unified.logRosenbluthReverse = safeLog(rosen);

    unified.logExtraForward = 0.0;
    unified.logExtraReverse = 0.0;

    auto eval = evaluate(unified, MoveType::DELETION);
    if (logRatioOut) {
        *logRatioOut = eval.logRatio;
    }
    return eval.probability;
}

// Calculate translation probability
double GCMCAcceptance::calculateTranslationProbability(
    double deltaE,
    double bias) {
    
    double beta = getBeta();
    return std::min(1.0, std::exp(-beta * deltaE) * bias);
}

// Calculate rotation probability
double GCMCAcceptance::calculateRotationProbability(
    double deltaE,
    double bias) {
    
    double beta = getBeta();
    return std::min(1.0, std::exp(-beta * deltaE) * bias);
}

// Calculate swap probability
double GCMCAcceptance::calculateSwapProbability(
    int type1, int type2,
    int n1, int n2,
    double deltaE,
    double bias) {
    
    double beta = getBeta();
    double activity1 = activities_[type1];
    double activity2 = activities_[type2];
    
    // Swap acceptance: min(1, (a2*n1)/(a1*(n2+1)) * exp(-beta*deltaE) * bias)
    double prefactor = (activity2 * n1) / (activity1 * (n2 + 1));
    double boltzmann = std::exp(-beta * deltaE);
    
    return std::min(1.0, prefactor * boltzmann * bias);
}

// Calculate volume change probability
double GCMCAcceptance::calculateVolumeChangeProbability(
    double oldVolume,
    double newVolume,
    int nMolecules,
    double deltaE) {
    
    double beta = getBeta();
    
    // NPT acceptance: min(1, (V_new/V_old)^N * exp(-beta*(deltaE + P*deltaV)))
    double volumeRatio = newVolume / oldVolume;
    double volumeTerm = std::pow(volumeRatio, nMolecules);
    double deltaV = newVolume - oldVolume;
    double pvTerm = pressure_ * deltaV;
    double boltzmann = std::exp(-beta * (deltaE + pvTerm));
    
    return std::min(1.0, volumeTerm * boltzmann);
}

// General acceptance calculation
double GCMCAcceptance::calculateAcceptance(
    CriterionType criterion,
    double deltaE,
    double bias,
    double additionalFactor) {
    
    double beta = getBeta();
    
    switch (criterion) {
        case CriterionType::METROPOLIS:
            return std::min(1.0, std::exp(-beta * deltaE) * bias * additionalFactor);
            
        case CriterionType::GRAND_CANONICAL:
            return std::min(1.0, std::exp(-beta * deltaE) * bias * additionalFactor);
            
        case CriterionType::ISOTHERMAL_ISOBARIC:
            return std::min(1.0, std::exp(-beta * deltaE) * bias * additionalFactor);
            
        case CriterionType::GIBBS_ENSEMBLE:
            return std::min(1.0, std::exp(-beta * deltaE) * bias * additionalFactor);
            
        case CriterionType::WANG_LANDAU:
            return calculateWangLandauAcceptance(1.0, additionalFactor);
            
        case CriterionType::TRANSITION_MATRIX:
            return std::min(1.0, bias * additionalFactor);
            
        default:
            return std::min(1.0, std::exp(-beta * deltaE) * bias * additionalFactor);
    }
}

// Accept or reject based on probability
bool GCMCAcceptance::acceptMove(double probability) {
    totalDecisions_++;
    double random = uniform_(rng_);
    bool accepted = (random < probability);
    if (accepted) acceptedMoves_++;
    return accepted;
}

// Accept using Metropolis criterion
bool GCMCAcceptance::acceptMetropolis(double deltaE, double bias) {
    double probability = calculateAcceptance(CriterionType::METROPOLIS, deltaE, bias);
    return acceptMove(probability);
}

// Calculate Gibbs ensemble acceptance
double GCMCAcceptance::calculateGibbsAcceptance(
    int boxFrom, int boxTo,
    int typeId,
    double deltaE) {
    
    // Suppress unused parameter warnings
    (void)boxFrom;
    (void)boxTo;
    (void)typeId;
    
    double beta = getBeta();
    // Simplified Gibbs acceptance
    return std::min(1.0, std::exp(-beta * deltaE));
}

// Calculate Wang-Landau acceptance
double GCMCAcceptance::calculateWangLandauAcceptance(
    double currentBias,
    double newBias) {
    
    // Wang-Landau: always accept if going to less visited state
    return std::min(1.0, currentBias / newBias);
}

// Calculate transition matrix acceptance
double GCMCAcceptance::calculateTransitionMatrixAcceptance(
    int oldState,
    int newState,
    const std::vector<std::vector<double>>& transitionMatrix) {
    
    if (oldState >= 0 && oldState < static_cast<int>(transitionMatrix.size()) &&
        newState >= 0 && newState < static_cast<int>(transitionMatrix[oldState].size())) {
        return transitionMatrix[oldState][newState];
    }
    return 0.0;
}

// Check detailed balance
bool GCMCAcceptance::checkDetailedBalance(
    double forwardProb,
    double reverseProb,
    double tolerance) {
    
    double ratio = forwardProb / (reverseProb + 1e-10);
    return std::abs(ratio - 1.0) < tolerance;
}

// Get ideal gas contribution
double GCMCAcceptance::getIdealGasContribution(int n, double V) const {
    if (n <= 0) return 0.0;
    // Stirling's approximation: ln(n!) ≈ n*ln(n) - n
    return n * std::log(V / n) + n;
}

// Get de Broglie wavelength
double GCMCAcceptance::getDeBroglieWavelength(double mass) const {
    // h = 6.626e-34 J·s, k = 1.381e-23 J/K
    const double h = 6.626e-34;
    const double k = 1.381e-23;
    const double pi = 3.14159265359;
    
    double lambda = h / std::sqrt(2.0 * pi * mass * k * temperature_);
    return lambda * 1e10; // Convert to Angstroms
}

double GCMCAcceptance::safeLog(double value) const {
    return std::log(std::max(value, 1e-30));
}

double GCMCAcceptance::safeExp(double logValue) const {
    if (logValue >= 0.0) {
        return 1.0;
    }
    const double lowerBound = -700.0;
    return std::exp(std::max(logValue, lowerBound));
}

gcmc::GrandCanonicalEvaluation GCMCAcceptance::evaluate(
    const gcmc::GrandCanonicalTerms& terms,
    MoveType moveType) {
    
    auto safeLogLocal = [](double value) {
        return std::log(std::max(value, 1e-30));
    };
    auto safeExpLocal = [](double logValue) {
        if (logValue >= 0.0) {
            return 1.0;
        }
        const double lowerBound = -700.0;
        return std::exp(std::max(logValue, lowerBound));
    };

    gcmc::GrandCanonicalEvaluation eval;
    
    if (moveType == MoveType::INSERTION) {
        int countAfter = terms.countAfter > 0 ? terms.countAfter : terms.countBefore + 1;
        if (countAfter <= 0) {
            eval.probability = 0.0;
            eval.logRatio = -std::numeric_limits<double>::infinity();
            return eval;
        }
        double logNplus1 = safeLogLocal(static_cast<double>(countAfter));
        double logRatio =
            (terms.logProposalReverse - terms.logProposalForward)
            - terms.beta * terms.deltaEnergy
            + terms.beta * terms.chemicalPotential
            + terms.logVolume
            + (terms.logCavityForward - terms.logCavityReverse)
            - logNplus1
            + (terms.logRosenbluthForward - terms.logRosenbluthReverse)
            - terms.logLambda3
            + (terms.logExtraForward - terms.logExtraReverse);
        eval.logRatio = logRatio;
        eval.probability = safeExpLocal(logRatio);
        return eval;
    }

    if (terms.countBefore <= 0) {
        eval.probability = 0.0;
        eval.logRatio = -std::numeric_limits<double>::infinity();
        return eval;
    }

	double logN = safeLogLocal(static_cast<double>(terms.countBefore));
	double logRatio =
	    (terms.logProposalReverse - terms.logProposalForward)
	    - terms.beta * terms.deltaEnergy
	    - terms.beta * terms.chemicalPotential
	    + logN
	    - terms.logVolume
	    + (terms.logCavityForward - terms.logCavityReverse)
        + (terms.logRosenbluthForward - terms.logRosenbluthReverse)
        + terms.logLambda3
        + (terms.logExtraForward - terms.logExtraReverse);
    eval.logRatio = logRatio;
    eval.probability = safeExpLocal(logRatio);
    return eval;
}

double GCMCAcceptance::getChemicalPotentialInternal(int typeId) const {
    auto muIt = chemicalPotentials_.find(typeId);
    if (muIt != chemicalPotentials_.end()) {
        return muIt->second;
    }
    auto actIt = activities_.find(typeId);
    if (actIt != activities_.end()) {
        double activity = std::max(actIt->second, 1e-30);
        double beta = getBeta();
        if (beta > 0.0) {
            return std::log(activity) / beta;
        }
    }
    return 0.0;
}

// ============================================================================
// GCMCSmartAcceptance Implementation
// ============================================================================

GCMCSmartAcceptance::GCMCSmartAcceptance()
    : GCMCAcceptance(),
      learningRate_(0.01) {
    // Initialize weights
    weights_.resize(10, 0.1);
}

// Learn from history
void GCMCSmartAcceptance::learnFromHistory(
    const std::vector<double>& deltaEs,
    const std::vector<bool>& accepted,
    const std::vector<double>& systemProperties) {
    
    // Simple gradient descent learning
    for (size_t i = 0; i < deltaEs.size(); ++i) {
        std::vector<double> features = extractFeatures(deltaEs[i], systemProperties);
        double prediction = predictAcceptance(deltaEs[i], systemProperties);
        double error = (accepted[i] ? 1.0 : 0.0) - prediction;
        
        // Update weights
        for (size_t j = 0; j < weights_.size() && j < features.size(); ++j) {
            weights_[j] += learningRate_ * error * features[j];
        }
    }
}

// Predict acceptance probability
double GCMCSmartAcceptance::predictAcceptance(
    double deltaE,
    const std::vector<double>& systemProperties) {
    
    std::vector<double> features = extractFeatures(deltaE, systemProperties);
    
    // Linear combination of features
    double score = 0.0;
    for (size_t i = 0; i < weights_.size() && i < features.size(); ++i) {
        score += weights_[i] * features[i];
    }
    
    // Sigmoid activation
    return 1.0 / (1.0 + std::exp(-score));
}

// Adjust criterion dynamically
void GCMCSmartAcceptance::adjustCriterion(double targetAcceptance) {
    double currentAcceptance = getAverageAcceptance();
    
    if (currentAcceptance < targetAcceptance - 0.1) {
        // Too low acceptance - increase temperature artificially
        temperature_ *= 1.05;
    } else if (currentAcceptance > targetAcceptance + 0.1) {
        // Too high acceptance - decrease temperature artificially
        temperature_ *= 0.95;
    }
}

// Extract features for machine learning
std::vector<double> GCMCSmartAcceptance::extractFeatures(
    double deltaE,
    const std::vector<double>& properties) {
    
    std::vector<double> features;
    
    // Energy-based features
    features.push_back(deltaE);
    features.push_back(deltaE * deltaE);
    features.push_back(std::exp(-getBeta() * deltaE));
    
    // System property features
    for (double prop : properties) {
        features.push_back(prop);
    }
    
    // Ensure we have at least 10 features
    while (features.size() < 10) {
        features.push_back(0.0);
    }
    
    return features;
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
