#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_ACTIVE_POOL_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_ACTIVE_POOL_HPP

#include <vector>
#include <queue>
#include <set>
#include <memory>
#include <utility>
#include "../common/MovementUtils.hpp"

// Forward declarations
namespace pygcmc {
namespace model {
namespace montecarlo {
    struct MCAtom;
    struct MCResidue;
    struct MCState;
}
}
}

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCAtom;
using model::montecarlo::MCResidue;
using model::montecarlo::MCState;

/**
 * Active Pool manages pre-allocated memory for atoms and residues
 * Uses active mask approach for efficient insertion/deletion
 */
class ActivePool {
public:
    // Residue metadata (since MCResidue is read-only)
    struct ResidueMetadata {
        int residueType = 0;
        bool active = false;
        int atomStartIndex = -1;
        int atomCount = 0;
        Vector3 centerOfMass;
        double insertionTime = 0.0;  // For analysis
    };
    
    // Statistics for monitoring
    struct Statistics {
        int totalInserts = 0;
        int totalDeletes = 0;
        int compactions = 0;
        int peakAtoms = 0;
        int peakResidues = 0;
        double averageFragmentation = 0.0;
        int batchOperations = 0;
    };
    
    // Constructor
    explicit ActivePool(int maxAtoms = 100000, int maxResidues = 30000);
    
    // Destructor
    ~ActivePool();
    
    // Core operations
    int insertMolecule(const std::vector<MCAtom>& atoms, int resType = 0);
    bool deleteResidue(int resIdx);
    int compact(bool force = false);
    
    // Batch operations for efficiency
    void queueInsert(const std::vector<MCAtom>& atoms, int resType = 0);
    void queueDelete(int resIdx);
    std::pair<int, int> flushBatch();  // Returns (inserted, deleted)
    
    // State synchronization
    void syncToState(MCState& state);
    void syncFromState(const MCState& state);
    
    // Query operations
    double getFragmentation() const;
    std::pair<int, int> getActiveCounts() const;  // (activeAtoms, activeResidues)
    std::vector<int> getActiveResidueIndices() const;
    bool isResidueActive(int resIdx) const;
    
    // Metadata access
    ResidueMetadata* getResidueMetadata(int resIdx);
    const ResidueMetadata* getResidueMetadata(int resIdx) const;
    
    // Statistics
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
    // Configuration
    void setFragmentationThreshold(double threshold) { fragmentationThreshold_ = threshold; }
    double getFragmentationThreshold() const { return fragmentationThreshold_; }
    
    // Capacity management
    int getMaxAtoms() const { return maxAtoms_; }
    int getMaxResidues() const { return maxResidues_; }
    bool canInsert(int atomCount) const;
    
private:
    // Pre-allocated arrays
    std::vector<MCAtom> atoms_;                  // Pre-allocated atom array
    std::vector<ResidueMetadata> residueInfo_;   // Residue metadata
    std::vector<bool> atomActive_;               // Active mask for atoms
    std::vector<bool> residueActive_;            // Active mask for residues
    
    // Free slot management
    std::set<int> freeResidueSlots_;            // Available residue indices
    std::queue<int> freeAtomRanges_;            // Available atom index ranges
    
    // Batch operation queues
    struct InsertOperation {
        std::vector<MCAtom> atoms;
        int resType;
    };
    std::queue<InsertOperation> insertQueue_;
    std::queue<int> deleteQueue_;
    
    // Capacity and thresholds
    int maxAtoms_;
    int maxResidues_;
    double fragmentationThreshold_;
    
    // Current state
    int activeAtomCount_;
    int activeResidueCount_;
    int nextFreeAtomIndex_;
    int nextFreeResidueIndex_;
    
    // Statistics
    mutable Statistics stats_;
    
    // Helper functions
    int findFreeResidueSlot();
    int findFreeAtomRange(int count);
    void updateFragmentationStats();
    void performCompaction();
    bool shouldCompact() const;
    
    // Copy atoms to pool
    void copyAtomsToPool(const std::vector<MCAtom>& atoms, int startIdx);
    
    // Update center of mass for a residue
    void updateCenterOfMass(int resIdx);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_ACTIVE_POOL_HPP