#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ENGINE_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ENGINE_HPP

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../reservoir/fragment_reservoir.hpp"
#include "../bias/CavityBias.hpp"
#include "../bias/ConfigBias.hpp"
#include "../../energy/EnergyModule.hpp"
#include "GCMCEnergyCallback.hpp"
#include "GCMCAcceptance.hpp"
#include "GCMCStatistics.hpp"
#include "../common/RegionConstraint.hpp"
#include <random>
#include <memory>
#include <unordered_map>
#include <string>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

using namespace model::montecarlo;
using namespace energy;

/**
 * @brief Core GCMC engine that handles all move operations
 * 
 * This class implements the core GCMC algorithm, managing:
 * - Fragment insertion and deletion
 * - Translation and rotation moves
 * - Energy calculations
 * - Acceptance/rejection decisions
 */
class GCMCEngine {
public:
    // Move result structure
    struct MoveResult {
        enum MoveType {
            INSERT,
            DELETE,
            TRANSLATE,
            ROTATE,
            SWAP,
            REGROWTH,
            CLUSTER
        };
        
        MoveType type;
        bool accepted;
        double energyBefore;
        double energyAfter;
        double deltaE;
        double bias;
        double acceptanceProbability;  // Added for observability
        int fragmentType;
        int residueIndex;
        Vector3 position;
        
        MoveResult() : type(INSERT), accepted(false), energyBefore(0), 
                      energyAfter(0), deltaE(0), bias(1.0), 
                      acceptanceProbability(0.0), fragmentType(-1), residueIndex(-1) {}
    };
    
    // Constructor
    GCMCEngine();
    ~GCMCEngine();
    
    // Initialize
    void initialize(MCState* state, FragmentReservoir* reservoir);
    
    // Set components
    void setCavityManager(CavityManager* cavityManager) {
        cavityManager_ = cavityManager;
    }
    void setConfigBiasManager(ConfigBiasManager* configBias) {
        configBias_ = configBias;
    }
    void setRegionConstraint(std::unique_ptr<RegionConstraint> constraint) {
        regionConstraint_ = std::move(constraint);
    }
    
    // Core move operations
    MoveResult attemptInsertion(int typeId);
    MoveResult attemptDeletion(int typeId);
    MoveResult attemptTranslation(int residueIdx);
    MoveResult attemptRotation(int residueIdx);
    MoveResult attemptSwap(int typeId1, int typeId2);
    MoveResult attemptRegrowth(int residueIdx);
    MoveResult attemptClusterMove(int residueIdx, double cutoff);
    
    // Fragment selection
    int selectRandomFragment();
    int selectRandomInstance(int typeId = -1);
    std::vector<int> selectCluster(int seedIdx, double cutoff);
    
    // Position and orientation generation
    Vector3 generateRandomPosition();
    Vector3 generateCavityPosition();
    Quaternion generateRandomOrientation();
    Vector3 generateTranslationVector(double maxDist);
    Quaternion generateRotationQuaternion(double maxAngle);
    
    // Energy calculations
    double calculateSystemEnergy();
    double calculateFragmentEnergy(int residueIdx);
    double calculateInteractionEnergy(int residueIdx);
    double calculatePairEnergy(int idx1, int idx2);
    
    // Bias calculations
    double calculateInsertionBias(const FragmentTemplate& tmpl, 
                                  const Vector3& position,
                                  const Quaternion& orientation);
    double calculateDeletionBias(int residueIdx);
    double calculateDeletionBiasAtPosition(const Vector3& position);
    
    // State synchronization
    void synchronizeStateWithReservoir(int instanceId, bool isInsertion);
    void updateAtomCoordinates(int residueIdx);  // Update positions without changing count
    double calculateRegrowthBias(int residueIdx);
    
    // Acceptance criteria
    bool acceptMove(double deltaE, double bias, double temperature);
    double calculateAcceptanceProbability(const MoveResult& result, 
                                         double temperature);
    
    // Residue position and orientation queries
    Vector3 getResiduePosition(int residueIdx);
    Quaternion getResidueOrientation(int residueIdx);
    
    // Configuration
    void setTemperature(double T) { temperature_ = T; }
    void setEnergyMethod(EnergyMethod method) { 
        energyMethod_ = method; 
        if (energyCallback_) {
            energyCallback_->setEnergyMethod(method);
        }
    }
    void setCutoff(double cutoff) { cutoff_ = cutoff; }
    void setSeed(unsigned int seed);
    void setAcceptanceCalculator(GCMCAcceptance* acceptCalc) { 
        acceptanceCalculator_ = acceptCalc;
        // Auto-seed acceptance if engine has been seeded
        if (acceptCalc && lastSeed_ != 0) {
            acceptCalc->setSeed(lastSeed_ + 1);
        }
    }
    
    // Set energy callback
    void setEnergyCallback(std::unique_ptr<GCMCEnergyCallback> callback) {
        energyCallback_ = std::move(callback);
    }

    // Get energy callback (for configuration)
    GCMCEnergyCallback* getEnergyCallback() {
        return energyCallback_.get();
    }

    // Set CBMC trials per fragment type
    void setCBMCTrialsPerType(const std::vector<int>& trials) {
        cbmcTrialsPerType_ = trials;
    }
    
    // Statistics
    int getTotalMoves() const { return totalMoves_; }
    int getAcceptedMoves() const { return acceptedMoves_; }
    double getAcceptanceRate() const {
        return totalMoves_ > 0 ? 
               static_cast<double>(acceptedMoves_) / totalMoves_ : 0.0;
    }
    
    // Dynamic configuration
    void setConfigValue(const std::string& key, double value);
    double getConfigValue(const std::string& key) const;
    
    // Statistics management
    void enableStatistics(bool enable) {
        collectStats_ = enable;
        // When enabling, initialize countdown to current interval to avoid
        // immediate sampling and eliminate per-move modulo overhead
        if (collectStats_) {
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    void setStatisticsInterval(int interval);
    GCMCStatistics& getStatistics() { return statistics_; }
    const GCMCStatistics& getStatistics() const { return statistics_; }
    
private:
    // State and components
    MCState* state_;
    FragmentReservoir* reservoir_;
    CavityManager* cavityManager_;
    ConfigBiasManager* configBias_;
    GCMCAcceptance* acceptanceCalculator_;
    std::unique_ptr<GCMCEnergyCallback> energyCallback_;
    std::unique_ptr<RegionConstraint> regionConstraint_;
    
    // Parameters
    double temperature_;
    double cutoff_;
    EnergyMethod energyMethod_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    
    // Statistics
    int totalMoves_;
    int acceptedMoves_;
    
    // Configuration parameters for moves
    double maxTranslationStep_;
    double maxRotationAngleRad_;
    bool useCavityBias_;
    bool useConfBias_ = false;
    std::vector<int> cbmcTrialsPerType_;  // CBMC trials per fragment type
    unsigned int lastSeed_ = 0;  // Store last seed for auto-seeding acceptance
    
    // Smart statistics collection
    GCMCStatistics statistics_;
    bool collectStats_ = false;
    int statsInterval_ = 1000;  // Default: sample every 1000 steps
    int statsCountdown_ = 1000; // Lightweight countdown to next sample
    
    // Dynamic configuration
    std::unordered_map<std::string, double> configMap_;
    
    // Performance optimization: cache probability storage flag
    mutable bool storeProbabilityCached_ = false;
    mutable bool storeProbabilityValue_ = false;
    
    // Helper methods
    void updateFragmentPosition(int residueIdx, const Vector3& newPos);
    void updateFragmentOrientation(int residueIdx, const Quaternion& newOrient);
    void applyPeriodicBoundary(Vector3& position);
    double minimumImageDistance(const Vector3& r1, const Vector3& r2);
    
    // Probability storage control
    bool shouldStoreProbability() const;

    // CBMC helper methods
    struct TrialConfiguration {
        Vector3 position;
        Quaternion orientation;
        double energy;
        double weight;
    };
    TrialConfiguration performCBMCInsertion(int typeId, int numTrials);
    double calculateCBMCBias(const std::vector<TrialConfiguration>& trials, int selectedIdx);
    
    // Energy caching (optional optimization)
    struct EnergyCache {
        bool valid;
        double totalEnergy;
        std::vector<double> fragmentEnergies;
        
        void invalidate() { valid = false; }
    };
    EnergyCache energyCache_;
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ENGINE_HPP
