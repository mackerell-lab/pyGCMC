#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_BIAS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_BIAS_HPP

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../bias/CavityBias.hpp"
#include "../bias/ConfigBias.hpp"
#include "../reservoir/fragment_reservoir.hpp"
#include <memory>
#include <vector>
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Use types from the included MCMain.hpp
using pygcmc::model::montecarlo::MCState;
using pygcmc::model::montecarlo::Vector3;
using pygcmc::model::montecarlo::Quaternion;
using pygcmc::platform::cpu::movement::FragmentTemplate;

/**
 * @brief Calculates various biases for GCMC moves
 * 
 * This class manages the calculation of different biasing techniques:
 * - Cavity bias for intelligent insertion positions
 * - Configurational bias (CBMC) for molecular conformations
 * - Orientational bias for rigid molecules
 * - Distance-based biases for specific interactions
 */
class GCMCBias {
public:
    // Bias result structure
    struct BiasResult {
        double totalBias;
        double cavityBias;
        double configBias;
        double orientBias;
        double distanceBias;
        
        // Detailed information
        Vector3 selectedPosition;
        Quaternion selectedOrientation;
        std::vector<Vector3> trialPositions;
        std::vector<double> trialWeights;
        int selectedTrial;
        
        BiasResult() : totalBias(1.0), cavityBias(1.0), configBias(1.0),
                      orientBias(1.0), distanceBias(1.0), selectedTrial(0) {}
        
        double getCombinedBias() const {
            return cavityBias * configBias * orientBias * distanceBias;
        }
    };
    
    // Constructor
    GCMCBias();
    ~GCMCBias();
    
    // Initialize with components
    void initialize(MCState* state);
    void setCavityManager(CavityManager* cavityManager);
    void setConfigBiasManager(ConfigBiasManager* configBias);
    
    // Calculate insertion bias
    BiasResult calculateInsertionBias(const FragmentTemplate& tmpl,
                                      int nTrials = 10);
    
    // Calculate deletion bias (reverse of insertion)
    BiasResult calculateDeletionBias(int residueIdx,
                                     const FragmentTemplate& tmpl,
                                     int nTrials = 10);
    
    // Calculate regrowth bias
    BiasResult calculateRegrowthBias(int residueIdx,
                                     const FragmentTemplate& tmpl,
                                     int nTrials = 10);
    
    // Individual bias components
    double calculateCavityBias(const Vector3& position);
    double calculateConfigBias(const FragmentTemplate& tmpl,
                              const Vector3& position,
                              const Quaternion& orientation,
                              int nTrials);
    double calculateOrientationalBias(const FragmentTemplate& tmpl,
                                      const Vector3& position,
                                      int nTrials);
    double calculateDistanceBias(const Vector3& position,
                                 const Vector3& target,
                                 double sigma);
    
    // Advanced biasing methods
    double calculatePreferentialSampling(const FragmentTemplate& tmpl,
                                         const std::vector<Vector3>& hotspots);
    double calculateUmbrellaSampling(double currentValue,
                                     double targetValue,
                                     double force);
    
    // Rosenbluth weight calculation
    double calculateRosenbluthWeight(const std::vector<double>& energies,
                                     double temperature);
    int selectRosenbluthTrial(const std::vector<double>& weights);
    
    // Configuration
    void setTemperature(double T) { temperature_ = T; }
    void enableCavityBias(bool enable) { useCavityBias_ = enable; }
    void enableConfigBias(bool enable) { useConfigBias_ = enable; }
    void enableOrientBias(bool enable) { useOrientBias_ = enable; }
    
    // Adaptive biasing
    void enableAdaptiveBiasing();
    void updateBiasParameters();
    
    // Statistics
    double getAverageBias() const { return averageBias_; }
    int getBiasCalculations() const { return biasCalculations_; }
    
private:
    // State and components
    MCState* state_;
    CavityManager* cavityManager_;
    ConfigBiasManager* configBias_;
    
    // Parameters
    double temperature_;
    bool useCavityBias_;
    bool useConfigBias_;
    bool useOrientBias_;
    bool useDistanceBias_;
    
    // Adaptive parameters
    bool adaptiveBiasing_;
    double biasPotential_;
    std::vector<double> biasHistory_;
    
    // Statistics
    double averageBias_;
    int biasCalculations_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;
    
    // Helper methods
    std::vector<Vector3> generateTrialPositions(int nTrials);
    std::vector<Quaternion> generateTrialOrientations(int nTrials);
    double evaluatePosition(const Vector3& position,
                           const FragmentTemplate& tmpl);
    double evaluateOrientation(const Quaternion& orientation,
                              const Vector3& position,
                              const FragmentTemplate& tmpl);
};

/**
 * @brief Wang-Landau bias calculator for flat histogram sampling
 */
class GCMCWangLandauBias : public GCMCBias {
public:
    GCMCWangLandauBias();
    
    // Initialize histogram
    void initializeHistogram(double minValue, double maxValue, int nBins);
    
    // Update histogram and bias
    void updateHistogram(double value);
    double getBias(double value) const;
    
    // Check convergence
    bool isConverged() const;
    
private:
    std::vector<double> histogram_;
    std::vector<double> biasFunction_;
    double minValue_, maxValue_;
    int nBins_;
    double modificationFactor_;
    double convergenceCriterion_;
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_BIAS_HPP