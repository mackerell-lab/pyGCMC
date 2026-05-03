#include "GCMCBias.hpp"
#include "../reservoir/FragmentReservoir.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Constructor
GCMCBias::GCMCBias()
    : state_(nullptr),
      cavityManager_(nullptr),
      configBias_(nullptr),
      temperature_(300.0),
      useCavityBias_(true),
      useConfigBias_(false),
      useOrientBias_(false),
      useDistanceBias_(false),
      adaptiveBiasing_(false),
      biasPotential_(0.0),
      averageBias_(1.0),
      biasCalculations_(0),
      rng_(0u),  // Deterministic initial value, will be set via setSeed()
      uniform_(0.0, 1.0),
      normal_(0.0, 1.0) {
}

// Destructor
GCMCBias::~GCMCBias() {
}

// Initialize with state
void GCMCBias::initialize(MCState* state) {
    state_ = state;
}

// Set cavity manager
void GCMCBias::setCavityManager(CavityManager* cavityManager) {
    cavityManager_ = cavityManager;
}

// Set config bias manager
void GCMCBias::setConfigBiasManager(ConfigBiasManager* configBias) {
    configBias_ = configBias;
}

// Calculate insertion bias
GCMCBias::BiasResult GCMCBias::calculateInsertionBias(
    const FragmentTemplate& tmpl, int nTrials) {

    BiasResult result;

    // Generate trial positions
    result.trialPositions = generateTrialPositions(nTrials);
    result.trialWeights.resize(nTrials);

    // Calculate weights for each trial
    double totalWeight = 0.0;
    for (int i = 0; i < nTrials; ++i) {
        double weight = 1.0;

        // Cavity bias
        if (useCavityBias_) {
            weight *= calculateCavityBias(result.trialPositions[i]);
        }

        // Evaluation energy (simplified)
        weight *= std::exp(-evaluatePosition(result.trialPositions[i], tmpl) / (8.314e-3 * temperature_));

        result.trialWeights[i] = weight;
        totalWeight += weight;
    }

    // Select trial based on weights
    double random = uniform_(rng_) * totalWeight;
    double cumWeight = 0.0;

    for (int i = 0; i < nTrials; ++i) {
        cumWeight += result.trialWeights[i];
        if (random <= cumWeight) {
            result.selectedTrial = i;
            result.selectedPosition = result.trialPositions[i];
            break;
        }
    }

    // Calculate bias factors
    result.cavityBias = useCavityBias_ ?
        calculateCavityBias(result.selectedPosition) : 1.0;

    // Generate and select orientation
    if (useOrientBias_) {
        result.orientBias = calculateOrientationalBias(tmpl, result.selectedPosition, nTrials);
    }

    // Config bias (if enabled)
    if (useConfigBias_) {
        result.configBias = calculateConfigBias(tmpl, result.selectedPosition,
                                               result.selectedOrientation, nTrials);
    }

    result.totalBias = result.getCombinedBias();

    // Update statistics
    biasCalculations_++;
    averageBias_ = (averageBias_ * (biasCalculations_ - 1) + result.totalBias) / biasCalculations_;

    return result;
}

// Calculate deletion bias
GCMCBias::BiasResult GCMCBias::calculateDeletionBias(
    int residueIdx, const FragmentTemplate& tmpl, int nTrials) {

    // Suppress unused parameter warnings
    (void)residueIdx;
    (void)tmpl;
    (void)nTrials;

    BiasResult result;

    // For deletion, we need to calculate the reverse bias
    // (what would be the bias to insert at this position)

    if (!state_) return result;

    // Get current position (simplified - would need actual residue position)
    Vector3 currentPos(0, 0, 0);  // Placeholder

    // Calculate what the insertion bias would have been
    result.cavityBias = useCavityBias_ ?
        1.0 / calculateCavityBias(currentPos) : 1.0;

    // Config bias reverse
    if (useConfigBias_) {
        result.configBias = 1.0;  // Simplified
    }

    result.totalBias = result.getCombinedBias();

    return result;
}

// Calculate regrowth bias
GCMCBias::BiasResult GCMCBias::calculateRegrowthBias(
    int residueIdx, const FragmentTemplate& tmpl, int nTrials) {

    // Regrowth = deletion + insertion
    BiasResult deletionBias = calculateDeletionBias(residueIdx, tmpl, nTrials);
    BiasResult insertionBias = calculateInsertionBias(tmpl, nTrials);

    BiasResult result;
    result.totalBias = deletionBias.totalBias * insertionBias.totalBias;
    result.cavityBias = deletionBias.cavityBias * insertionBias.cavityBias;
    result.configBias = deletionBias.configBias * insertionBias.configBias;
    result.orientBias = deletionBias.orientBias * insertionBias.orientBias;
    result.selectedPosition = insertionBias.selectedPosition;
    result.selectedOrientation = insertionBias.selectedOrientation;

    return result;
}

// Calculate cavity bias
double GCMCBias::calculateCavityBias(const Vector3& position) {
    if (!cavityManager_) return 1.0;

    // Convert to movement::Vector3
    movement::Vector3 pos(position.x, position.y, position.z);

    // Check if position is in a cavity
    return cavityManager_->getCavityScore(pos);
}

// Calculate config bias
double GCMCBias::calculateConfigBias(
    const FragmentTemplate& tmpl,
    const Vector3& position,
    const Quaternion& orientation,
    int nTrials) {

    if (!configBias_ || !state_) return 1.0;

    // Generate trial configurations and calculate energies
    std::vector<double> trialEnergies;
    trialEnergies.reserve(nTrials);

    for (int i = 0; i < nTrials; ++i) {
        // Generate trial configuration by rotating the template
        double angle = uniform_(rng_) * 2 * M_PI;
        Vector3 axis(normal_(rng_), normal_(rng_), normal_(rng_));
        double axisNorm = axis.norm();
        if (axisNorm > 0) {
            axis = axis * (1.0 / axisNorm);  // Normalize manually
        }

        Quaternion trialOrientation = orientation;
        if (i > 0) {  // Keep first trial as original
            // Apply small random rotation
            double halfAngle = angle / 2;
            double s = std::sin(halfAngle);
            Quaternion rotation(std::cos(halfAngle), s * axis.x, s * axis.y, s * axis.z);
            // Quaternion multiplication: q1 * q2
            // (w1, x1, y1, z1) * (w2, x2, y2, z2) =
            // (w1*w2 - x1*x2 - y1*y2 - z1*z2,
            //  w1*x2 + x1*w2 + y1*z2 - z1*y2,
            //  w1*y2 - x1*z2 + y1*w2 + z1*x2,
            //  w1*z2 + x1*y2 - y1*x2 + z1*w2)
            trialOrientation.w = rotation.w * orientation.w - rotation.x * orientation.x - rotation.y * orientation.y - rotation.z * orientation.z;
            trialOrientation.x = rotation.w * orientation.x + rotation.x * orientation.w + rotation.y * orientation.z - rotation.z * orientation.y;
            trialOrientation.y = rotation.w * orientation.y - rotation.x * orientation.z + rotation.y * orientation.w + rotation.z * orientation.x;
            trialOrientation.z = rotation.w * orientation.z + rotation.x * orientation.y - rotation.y * orientation.x + rotation.z * orientation.w;
            trialOrientation.normalize();
        }

        // Calculate energy for this configuration
        // This would require temporarily placing the fragment and calculating energy
        // For now, use simplified energy based on orientation
        double energy = evaluateOrientation(trialOrientation, position, tmpl);
        trialEnergies.push_back(energy);
    }

    // Calculate and return Rosenbluth weight
    return calculateRosenbluthWeight(trialEnergies, temperature_);
}

// Calculate orientational bias
double GCMCBias::calculateOrientationalBias(
    const FragmentTemplate& tmpl,
    const Vector3& position,
    int nTrials) {

    std::vector<Quaternion> orientations = generateTrialOrientations(nTrials);
    std::vector<double> weights(nTrials);

    for (int i = 0; i < nTrials; ++i) {
        weights[i] = evaluateOrientation(orientations[i], position, tmpl);
    }

    return calculateRosenbluthWeight(weights, temperature_);
}

// Calculate distance bias
double GCMCBias::calculateDistanceBias(
    const Vector3& position,
    const Vector3& target,
    double sigma) {

    double dist2 = (position - target).norm2();
    return std::exp(-dist2 / (2.0 * sigma * sigma));
}

// Calculate preferential sampling
double GCMCBias::calculatePreferentialSampling(
    const FragmentTemplate& tmpl,
    const std::vector<Vector3>& hotspots) {

    // Suppress unused parameter warnings
    (void)tmpl;

    if (hotspots.empty()) return 1.0;

    // Find nearest hotspot
    double minDist2 = 1e10;
    for (size_t i = 0; i < hotspots.size(); ++i) {
        double dist2 = 0.0;  // Would calculate distance to template position
        minDist2 = std::min(minDist2, dist2);
    }

    return std::exp(-minDist2 / 100.0);  // Gaussian weight
}

// Calculate umbrella sampling bias
double GCMCBias::calculateUmbrellaSampling(
    double currentValue,
    double targetValue,
    double force) {

    double delta = currentValue - targetValue;
    return std::exp(-0.5 * force * delta * delta);
}

// Calculate Rosenbluth weight
double GCMCBias::calculateRosenbluthWeight(
    const std::vector<double>& energies,
    double temperature) {

    double beta = 1.0 / (8.314e-3 * temperature);
    double weight = 0.0;

    for (double energy : energies) {
        weight += std::exp(-beta * energy);
    }

    // Rosenbluth weight is the sum of Boltzmann factors, not the average
    return weight;
}

// Select Rosenbluth trial
int GCMCBias::selectRosenbluthTrial(const std::vector<double>& weights) {
    double totalWeight = std::accumulate(weights.begin(), weights.end(), 0.0);
    double random = uniform_(rng_) * totalWeight;
    double cumWeight = 0.0;

    for (size_t i = 0; i < weights.size(); ++i) {
        cumWeight += weights[i];
        if (random <= cumWeight) {
            return static_cast<int>(i);
        }
    }

    return static_cast<int>(weights.size() - 1);
}

// Enable adaptive biasing
void GCMCBias::enableAdaptiveBiasing() {
    adaptiveBiasing_ = true;
    biasHistory_.clear();
}

// Update bias parameters
void GCMCBias::updateBiasParameters() {
    if (!adaptiveBiasing_) return;

    // Simple adaptive scheme
    if (biasHistory_.size() > 100) {
        double avgRecent = std::accumulate(
            biasHistory_.end() - 50, biasHistory_.end(), 0.0) / 50.0;
        double avgOld = std::accumulate(
            biasHistory_.begin(), biasHistory_.begin() + 50, 0.0) / 50.0;

        // Adjust parameters based on trend
        if (avgRecent < avgOld * 0.8) {
            // Bias is decreasing - may need adjustment
            biasPotential_ *= 1.1;
        } else if (avgRecent > avgOld * 1.2) {
            // Bias is increasing
            biasPotential_ *= 0.9;
        }
    }
}

// Generate trial positions
std::vector<Vector3> GCMCBias::generateTrialPositions(int nTrials) {
    std::vector<Vector3> positions;

    if (!state_) {
        // Random positions in box
        for (int i = 0; i < nTrials; ++i) {
            Vector3 pos(
                uniform_(rng_) * 100.0 - 50.0,
                uniform_(rng_) * 100.0 - 50.0,
                uniform_(rng_) * 100.0 - 50.0
            );
            positions.push_back(pos);
        }
    } else {
        // Use box dimensions from state
        for (int i = 0; i < nTrials; ++i) {
            Vector3 pos(
                uniform_(rng_) * state_->periodicBox[0] - state_->periodicBox[0]/2,
                uniform_(rng_) * state_->periodicBox[1] - state_->periodicBox[1]/2,
                uniform_(rng_) * state_->periodicBox[2] - state_->periodicBox[2]/2
            );
            positions.push_back(pos);
        }
    }

    return positions;
}

// Generate trial orientations
std::vector<Quaternion> GCMCBias::generateTrialOrientations(int nTrials) {
    std::vector<Quaternion> orientations;

    for (int i = 0; i < nTrials; ++i) {
        // Random quaternion
        double u1 = uniform_(rng_);
        double u2 = uniform_(rng_);
        double u3 = uniform_(rng_);

        Quaternion q(
            std::sqrt(1 - u1) * std::sin(2 * M_PI * u2),
            std::sqrt(1 - u1) * std::cos(2 * M_PI * u2),
            std::sqrt(u1) * std::sin(2 * M_PI * u3),
            std::sqrt(u1) * std::cos(2 * M_PI * u3)
        );
        q.normalize();
        orientations.push_back(q);
    }

    return orientations;
}

// Evaluate position
double GCMCBias::evaluatePosition(const Vector3& position, const FragmentTemplate& tmpl) {
    // Suppress unused parameter warning
    (void)tmpl;

    // Simplified energy evaluation
    // In practice, would calculate interaction energy with system
    double energy = 0.0;

    // Simple soft-core potential to avoid hard overlaps
    double r = position.norm();
    if (r < 0.1) r = 0.1;  // Avoid singularity

    // Soft repulsive potential at origin (in kJ/mol)
    if (r < 1.0) {  // Within 1 nm of origin
        energy = 10.0 * (1.0 - r);  // Linear repulsion
    }

    return energy;
}

// Evaluate orientation
double GCMCBias::evaluateOrientation(
    const Quaternion& orientation,
    const Vector3& position,
    const FragmentTemplate& tmpl) {

    // Suppress unused parameter warnings
    (void)orientation;
    (void)position;
    (void)tmpl;

    // Simplified orientation evaluation
    // In practice, would calculate orientation-dependent interactions
    return 1.0;
}

// ============================================================================
// GCMCWangLandauBias Implementation
// ============================================================================

GCMCWangLandauBias::GCMCWangLandauBias()
    : GCMCBias(),
      minValue_(0.0),
      maxValue_(1.0),
      nBins_(100),
      modificationFactor_(1.0),
      convergenceCriterion_(0.8) {
}

// Initialize histogram
void GCMCWangLandauBias::initializeHistogram(double minValue, double maxValue, int nBins) {
    minValue_ = minValue;
    maxValue_ = maxValue;
    nBins_ = nBins;

    histogram_.resize(nBins_, 0.0);
    biasFunction_.resize(nBins_, 0.0);
}

// Update histogram
void GCMCWangLandauBias::updateHistogram(double value) {
    if (value < minValue_ || value > maxValue_) return;

    int bin = static_cast<int>((value - minValue_) / (maxValue_ - minValue_) * nBins_);
    if (bin >= 0 && bin < nBins_) {
        histogram_[bin] += 1.0;
        biasFunction_[bin] += modificationFactor_;
    }
}

// Get bias
double GCMCWangLandauBias::getBias(double value) const {
    if (value < minValue_ || value > maxValue_) return 1.0;

    int bin = static_cast<int>((value - minValue_) / (maxValue_ - minValue_) * nBins_);
    if (bin >= 0 && bin < nBins_) {
        return std::exp(biasFunction_[bin]);
    }

    return 1.0;
}

// Check convergence
bool GCMCWangLandauBias::isConverged() const {
    if (histogram_.empty()) return false;

    double minCount = *std::min_element(histogram_.begin(), histogram_.end());
    double maxCount = *std::max_element(histogram_.begin(), histogram_.end());

    if (maxCount > 0) {
        return (minCount / maxCount) > convergenceCriterion_;
    }

    return false;
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
