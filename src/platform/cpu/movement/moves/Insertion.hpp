#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_INSERTION_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_INSERTION_HPP

#include <memory>
#include <vector>
#include "../common/MovementInterface.hpp"
#include "../common/MovementParams.hpp"
#include "../common/MovementResult.hpp"
#include "../common/MovementUtils.hpp"

// Forward declarations for MCAtom
namespace pygcmc {
namespace model {
namespace montecarlo {
    struct MCAtom;
}
}
}

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCAtom;

// Forward declarations
class ActivePool;
class CavityManager;
class CavityBiasInsertion;
class EnergyInterface;

/**
 * Insertion move implementation
 * Handles molecule insertion with optional cavity bias
 */
class InsertionMove : public MovementInterface {
public:
    // Constructor
    InsertionMove(ActivePool* activePool, 
                  CavityManager* cavityManager,
                  EnergyInterface* energyCalc);
    
    // Destructor
    virtual ~InsertionMove();
    
    // Perform insertion attempt
    virtual MovementResult attemptInsertion(MCState& state, const MovementParams& params) override;
    virtual MovementResult attemptDeletion(MCState& /*state*/, const MovementParams& /*params*/) override { 
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "delete");
    }
    virtual MovementResult attemptTranslation(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "translate");
    }
    virtual MovementResult attemptRotation(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "rotate");
    }
    
    // Specific insertion methods
    MovementResult performSimpleInsertion(MCState& state, const MovementParams& params, int moleculeType = 0);
    MovementResult performCavityBiasInsertion(MCState& state, const MovementParams& params, int moleculeType = 0);
    
    // Configuration
    void setMoleculeType(int type) { moleculeType_ = type; }
    int getMoleculeType() const { return moleculeType_; }
    
    // Statistics
    struct Statistics {
        int totalAttempts = 0;
        int cavityInsertions = 0;
        int randomInsertions = 0;
        int acceptedInsertions = 0;
        double averageEnergyChange = 0.0;
        double averageCavityBias = 0.0;
        double acceptanceRate() const {
            return totalAttempts > 0 ? static_cast<double>(acceptedInsertions) / totalAttempts : 0.0;
        }
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
protected:
    // Create a molecule at a given position
    std::vector<MCAtom> createMolecule(int moleculeType, const Vector3& position);
    
    // Calculate insertion acceptance probability
    double calculateInsertionProbability(
        int n,
        double deltaE,
        const MovementParams& params,
        double cavityBiasFactor
    );
    
    // Select insertion position
    Vector3 selectInsertionPosition(MCState& state, const MovementParams& params, double& cavityBias);
    
private:
    ActivePool* activePool_;
    CavityManager* cavityManager_;
    EnergyInterface* energyCalc_;
    std::unique_ptr<CavityBiasInsertion> cavityBiasInsertion_;
    
    int moleculeType_;
    Statistics stats_;
    
    // Update statistics
    void updateStatistics(bool accepted, bool usedCavity, double energyChange, double cavityBias);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_INSERTION_HPP