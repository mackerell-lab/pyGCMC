#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_ROTATION_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_ROTATION_HPP

#include "../common/MovementInterface.hpp"
#include "../common/MovementParams.hpp"
#include "../common/MovementResult.hpp"
#include "../common/MovementUtils.hpp"
#include <vector>
#include <memory>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Forward declarations
class ActivePool;
class ConfigBiasManager;
class EnergyInterface;

/**
 * Rotation move implementation
 * Handles molecule rotation with optional configurational bias
 */
class RotationMove : public MovementInterface {
public:
    // Constructor
    RotationMove(ActivePool* activePool, 
                 ConfigBiasManager* configBiasManager,
                 EnergyInterface* energyCalc);
    
    // Destructor
    virtual ~RotationMove();
    
    // Perform rotation attempt
    virtual MovementResult attemptRotation(MCState& state, const MovementParams& params) override;
    virtual MovementResult attemptInsertion(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "insert");
    }
    virtual MovementResult attemptDeletion(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "delete");
    }
    virtual MovementResult attemptTranslation(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "translate");
    }
    
    // Specific rotation methods
    MovementResult performSimpleRotation(MCState& state, const MovementParams& params, int residueIndex = -1);
    MovementResult performConfigBiasRotation(MCState& state, const MovementParams& params, int residueIndex = -1);
    
    // Statistics
    struct Statistics {
        int totalAttempts = 0;
        int simpleRotations = 0;
        int configBiasRotations = 0;
        int acceptedRotations = 0;
        int rejectedEmpty = 0;  // Rejected because no molecules to rotate
        double averageAngle = 0.0;
        double averageEnergyChange = 0.0;
        double averageBiasFactor = 0.0;
        double acceptanceRate() const {
            return totalAttempts > 0 ? static_cast<double>(acceptedRotations) / totalAttempts : 0.0;
        }
        double configBiasAcceptanceRate() const {
            return configBiasRotations > 0 ? 
                static_cast<double>(acceptedRotations - simpleRotations + acceptedRotations) / configBiasRotations : 0.0;
        }
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
protected:
    // Configuration for rotation
    struct RotationConfig {
        Quaternion quaternion;
        Vector3 center;  // Center of rotation
        std::vector<Vector3> originalPositions;
        double originalEnergy;
    };
    
    // Select a residue for rotation
    int selectResidueForRotation(const MCState& state);
    
    // Generate random rotation
    Quaternion generateRandomRotation(double maxAngle = M_PI);
    
    // Save configuration
    RotationConfig saveConfiguration(const MCState& state, int residueIndex);
    
    // Restore configuration
    void restoreConfiguration(MCState& state, int residueIndex, const RotationConfig& config);
    
    // Rotate a residue
    void rotateResidue(MCState& state, int residueIndex, const Quaternion& quaternion);
    
    // Calculate center of mass for rotation
    Vector3 calculateCenterOfMass(const MCState& state, int residueIndex);
    
    // Configurational bias rotation implementation
    struct ConfigBiasRotationResult {
        bool accepted;
        double energyChange;
        double biasFactor;
        int selectedConfigIndex;
        int totalConfigs;
    };
    
    ConfigBiasRotationResult performConfigBiasRotationInternal(
        MCState& state, 
        int residueIndex,
        const MovementParams& params
    );
    
private:
    ActivePool* activePool_;
    ConfigBiasManager* configBiasManager_;
    EnergyInterface* energyCalc_;
    Statistics stats_;
    
    // Update statistics
    void updateStatistics(bool accepted, bool usedConfigBias, 
                         double angle, double energyChange, double biasFactor);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_ROTATION_HPP