#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_DELETION_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_DELETION_HPP

#include "../common/MovementInterface.hpp"
#include "../common/MovementParams.hpp"
#include "../common/MovementResult.hpp"
#include "../common/MovementUtils.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Forward declarations
class ActivePool;
class EnergyInterface;

/**
 * Deletion move implementation
 * Handles molecule deletion from the system
 */
class DeletionMove : public MovementInterface {
public:
    // Constructor
    DeletionMove(ActivePool* activePool, EnergyInterface* energyCalc);
    
    // Destructor
    virtual ~DeletionMove();
    
    // Perform deletion attempt
    virtual MovementResult attemptDeletion(MCState& state, const MovementParams& params) override;
    virtual MovementResult attemptInsertion(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "insert");
    }
    virtual MovementResult attemptTranslation(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "translate");
    }
    virtual MovementResult attemptRotation(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "rotate");
    }
    
    // Specific deletion method
    MovementResult performDeletion(MCState& state, const MovementParams& params, int residueIndex = -1);
    
    // Statistics
    struct Statistics {
        int totalAttempts = 0;
        int acceptedDeletions = 0;
        int rejectedEmpty = 0;  // Rejected because no molecules to delete
        double averageEnergyChange = 0.0;
        double acceptanceRate() const {
            return totalAttempts > 0 ? static_cast<double>(acceptedDeletions) / totalAttempts : 0.0;
        }
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
protected:
    // Select a residue for deletion
    int selectResidueForDeletion(const MCState& state);
    
    // Calculate deletion acceptance probability
    double calculateDeletionProbability(
        int n,
        double deltaE,
        const MovementParams& params
    );
    
    // Calculate energy of a single residue
    double calculateResidueEnergy(const MCState& state, int residueIndex);
    
private:
    ActivePool* activePool_;
    EnergyInterface* energyCalc_;
    Statistics stats_;
    
    // Update statistics
    void updateStatistics(bool accepted, double energyChange);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_DELETION_HPP