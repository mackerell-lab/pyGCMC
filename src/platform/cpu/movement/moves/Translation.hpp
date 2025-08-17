#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_TRANSLATION_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_TRANSLATION_HPP

#include "../common/MovementInterface.hpp"
#include "../common/MovementParams.hpp"
#include "../common/MovementResult.hpp"
#include "../common/MovementUtils.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Forward declarations
class ActivePool;
class EnergyInterface;

/**
 * Translation move implementation
 * Handles molecule translation (displacement) in the system
 */
class TranslationMove : public MovementInterface {
public:
    // Constructor
    TranslationMove(ActivePool* activePool, EnergyInterface* energyCalc);
    
    // Destructor
    virtual ~TranslationMove();
    
    // Perform translation attempt
    virtual MovementResult attemptTranslation(MCState& state, const MovementParams& params) override;
    virtual MovementResult attemptInsertion(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "insert");
    }
    virtual MovementResult attemptDeletion(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "delete");
    }
    virtual MovementResult attemptRotation(MCState& /*state*/, const MovementParams& /*params*/) override {
        // Not implemented in this class
        return MovementResult(false, 0.0, 0.0, "rotate");
    }
    
    // Specific translation method
    MovementResult performTranslation(MCState& state, const MovementParams& params, int residueIndex = -1);
    
    // Batch translation for efficiency
    std::vector<MovementResult> performBatchTranslations(
        MCState& state, 
        const MovementParams& params, 
        int numAttempts
    );
    
    // Statistics
    struct Statistics {
        int totalAttempts = 0;
        int acceptedTranslations = 0;
        int rejectedEmpty = 0;  // Rejected because no molecules to translate
        double averageDisplacement = 0.0;
        double averageEnergyChange = 0.0;
        double acceptanceRate() const {
            return totalAttempts > 0 ? static_cast<double>(acceptedTranslations) / totalAttempts : 0.0;
        }
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
protected:
    // Select a residue for translation
    int selectResidueForTranslation(const MCState& state);
    
    // Generate random displacement
    Vector3 generateDisplacement(double maxTranslation);
    
    // Save atom positions for a residue
    std::vector<Vector3> saveAtomPositions(const MCState& state, int residueIndex);
    
    // Restore atom positions for a residue
    void restoreAtomPositions(MCState& state, int residueIndex, const std::vector<Vector3>& positions);
    
    // Translate a residue
    void translateResidue(MCState& state, int residueIndex, const Vector3& displacement);
    
    // Calculate energy before and after translation
    std::pair<double, double> calculateEnergyChange(
        MCState& state, 
        int residueIndex, 
        const Vector3& displacement
    );
    
private:
    ActivePool* activePool_;
    EnergyInterface* energyCalc_;
    Statistics stats_;
    
    // Update statistics
    void updateStatistics(bool accepted, double displacement, double energyChange);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_TRANSLATION_HPP