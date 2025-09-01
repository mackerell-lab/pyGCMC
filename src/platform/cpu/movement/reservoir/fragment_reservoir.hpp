// Fragment Reservoir System for PyGCMC
// Based on existing ActivePool and active/inactive mechanism
#pragma once

#include <vector>
#include <map>
#include <set>
#include <queue>
#include <memory>
#include <string>
#include <algorithm>
#include <cmath>
// Include MovementUtils for Vector3 and Quaternion first
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCStructures.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Forward declarations
class ActivePool;

// Import types from montecarlo namespace
using model::montecarlo::MCAtom;
using model::montecarlo::MCResidue;

// ============================================================================
// Fragment Template - Definition of insertable molecule types
// ============================================================================

struct FragmentTemplate {
    // Basic properties
    std::string name;                      // Fragment name (e.g., "WAT", "NA+")
    int typeId;                            // Type ID for force field
    std::vector<MCAtom> atoms;            // Atom template
    double molecularWeight;                // Molecular weight (g/mol)
    double radius;                         // Effective radius (Å)
    
    // Thermodynamic properties
    double chemicalPotential;             // Chemical potential μ (kJ/mol)
    double activity;                       // Activity z = exp(β*μ)
    double concentration;                  // Target concentration (M)
    
    // Topology information
    struct Bond {
        int atom1, atom2;
        double length;
    };
    struct Angle {
        int atom1, atom2, atom3;
        double angle;  // in radians
    };
    struct Dihedral {
        int atom1, atom2, atom3, atom4;
        double angle;  // in radians
    };
    
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    
    // Insertion preferences
    bool useCavityBias = true;            // Use cavity-biased insertion
    bool useConfigBias = false;           // Use configurational bias (CBMC)
    int configBiasTrials = 10;            // Number of trial configurations
    
    // Constraints
    bool isRigid = true;                  // Rigid molecule (no internal DOF)
    bool allowRotation = true;            // Allow rotation during MC
    bool allowTranslation = true;         // Allow translation during MC
    
    // Statistics (mutable for const access)
    mutable int totalInsertions = 0;
    mutable int successfulInsertions = 0;
    mutable double averageLifetime = 0.0;
    
    // Helper functions
    double getInsertionProbability() const {
        return totalInsertions > 0 ? 
            static_cast<double>(successfulInsertions) / totalInsertions : 0.0;
    }
    
    void updateActivity(double beta) {
        activity = std::exp(beta * chemicalPotential);
    }
    
    void calculateActivity(double temperature) {
        double beta = 1.0 / (8.314e-3 * temperature);  // kJ/mol/K
        updateActivity(beta);
    }
};

// ============================================================================
// Fragment Instance - Actual molecules in the system
// ============================================================================

struct FragmentInstance {
    // Identity
    int templateId;                        // Template ID
    int instanceId;                        // Unique instance ID
    int residueIndex;                      // Index in ActivePool
    
    // State
    bool isActive = true;                  // Active in simulation
    bool isGhost = false;                  // Ghost state (deleted but not recycled)
    bool isFixed = false;                  // Fixed position (no movement)
    
    // Position and orientation
    Vector3 centerOfMass;                  // Center of mass position
    Vector3 position;                       // Position (alias for centerOfMass)
    Quaternion orientation;                // Orientation quaternion
    Vector3 velocity;                      // Velocity (for future MD)
    
    // Energy
    double energy_vdw = 0.0;              // van der Waals energy
    double energy_elec = 0.0;             // Electrostatic energy  
    double energy_total = 0.0;            // Total energy
    double lastEnergyUpdate = 0.0;        // MC step of last energy update
    
    // Neighbor lists (cached)
    std::vector<int> proteinNeighbors;    // Neighboring protein atoms
    std::vector<int> fragmentNeighbors;   // Neighboring fragment instances
    double neighborListCutoff = 12.0;     // Cutoff for neighbor list
    double lastNeighborUpdate = 0.0;      // MC step of last update
    
    // History
    double insertionTime;                  // MC step when inserted
    double lastMoveTime;                   // MC step of last accepted move
    int moveAttempts = 0;                 // Total move attempts
    int acceptedMoves = 0;                // Accepted moves
    
    // Helper functions
    double getAcceptanceRate() const {
        return moveAttempts > 0 ? 
            static_cast<double>(acceptedMoves) / moveAttempts : 0.0;
    }
    
    double getLifetime(double currentStep) const {
        return currentStep - insertionTime;
    }
    
    bool needsNeighborUpdate(double currentStep, double updateFrequency = 100) const {
        return (currentStep - lastNeighborUpdate) > updateFrequency;
    }
};

// ============================================================================
// Fragment Reservoir - Main management class
// ============================================================================

class FragmentReservoir {
public:
    // Configuration
    struct Config {
        int maxInstances = 10000;          // Maximum fragment instances
        int maxGhosts = 100;               // Maximum ghost fragments to keep
        double ghostRecycleRatio = 0.5;    // Prefer ghost recycling when > 50% available
        bool autoCompact = true;           // Automatic memory compaction
        double compactThreshold = 0.3;     // Compact when fragmentation > 30%
        bool trackStatistics = true;       // Enable detailed statistics
        int statisticsWindow = 1000;       // Window for moving averages
    };
    
    // Statistics
    struct Statistics {
        // Per-type statistics
        std::map<int, int> activeCountByType;
        std::map<int, int> ghostCountByType;
        std::map<int, double> averageLifetimeByType;
        std::map<int, double> acceptanceRateByType;
        
        // Global statistics
        int totalInsertions = 0;
        int totalDeletions = 0;
        int ghostRecycles = 0;
        int memoryCompactions = 0;
        double averageGhostLifetime = 0.0;
        double peakActiveCount = 0;
        double averageActiveCount = 0.0;
        
        // Performance metrics
        double insertionTimeMs = 0.0;
        double deletionTimeMs = 0.0;
        double queryTimeMs = 0.0;
        int cacheHits = 0;
        int cacheMisses = 0;
        
        void print() const;
        void reset();
    };
    
public:
    // Constructor/Destructor
    explicit FragmentReservoir(std::shared_ptr<ActivePool> pool = nullptr);
    FragmentReservoir(const Config& config, std::shared_ptr<ActivePool> pool = nullptr);
    ~FragmentReservoir();
    
    // === Template Management ===
    
    // Add a new fragment template
    int addTemplate(const FragmentTemplate& tmpl);
    
    // Load template from file
    int loadTemplate(const std::string& filename, const std::string& name,
                    double chemicalPotential = 0.0);
    
    // Get template by ID or name
    FragmentTemplate* getTemplate(int templateId);
    const FragmentTemplate* getTemplate(int templateId) const;
    FragmentTemplate* getTemplate(const std::string& name);
    const FragmentTemplate* getTemplate(const std::string& name) const;
    int getTemplateCount() const { return templates_.size(); }
    
    // Update template properties
    void updateChemicalPotential(int templateId, double mu);
    void updateActivity(int templateId, double beta);
    void updateAllActivities(double beta);
    
    // === Instance Management ===
    
    // Create a new fragment instance
    int createInstance(int templateId, 
                      const Vector3& position,
                      const Quaternion& orientation = Quaternion());
    
    // Create with configurational bias
    int createInstanceCBMC(int templateId,
                          const std::vector<Vector3>& trialPositions,
                          const std::vector<double>& trialEnergies);
    
    // Delete an instance (convert to ghost)
    bool deleteInstance(int instanceId);
    
    // Restore a deleted (ghost) instance back to active state
    // Returns true if successfully restored, false otherwise
    bool restoreInstance(int instanceId, 
                        const Vector3& position,
                        const Quaternion& orientation);
    
    // Permanently remove an instance
    bool purgeInstance(int instanceId);
    
    // Batch operations
    std::vector<int> createMultipleInstances(int templateId, 
                                            const std::vector<Vector3>& positions);
    int deleteMultipleInstances(const std::vector<int>& instanceIds);
    
    // === Ghost Management ===
    
    // Recycle a ghost fragment
    int recycleGhost(int templateId);
    
    // Purge old ghosts
    int purgeGhosts(int maxToKeep = -1);
    
    // Get ghost statistics
    int getGhostCount(int templateId = -1) const;
    std::vector<int> getGhostIndices() const;
    
    // === Query Operations ===
    
    // Get active instances
    std::vector<int> getActiveInstances(int templateId = -1) const;
    int getActiveCount(int templateId = -1) const;
    
    // Get instance
    FragmentInstance* getInstance(int instanceId);
    const FragmentInstance* getInstance(int instanceId) const;
    FragmentInstance* getInstanceByResidueIndex(int residueIdx);
    
    // Find instances in region
    std::vector<int> findInstancesInSphere(const Vector3& center, double radius) const;
    std::vector<int> findInstancesInBox(const Vector3& min, const Vector3& max) const;
    
    // Get instances by property
    std::vector<int> getInstancesByEnergy(double minE, double maxE) const;
    std::vector<int> getInstancesByLifetime(double minTime, double maxTime, 
                                           double currentStep) const;
    
    // === Energy Management ===
    
    // Get energy (inline implementations at bottom)
    double getTotalEnergy(int instanceId) const;
    double getVdwEnergy(int instanceId) const;
    double getElecEnergy(int instanceId) const;
    
    // NOTE: Additional methods for neighbor management, movement tracking, 
    // and synchronization are not implemented in this stub version
    
    // Minimal methods required by GCMCEngine
    void updatePosition(int instanceId, const Vector3& newPos);
    void updateOrientation(int instanceId, const Quaternion& newOrient);
    
    // === Statistics ===
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    void printStatistics() const;
    
    // === Configuration ===
    
    const Config& getConfig() const { return config_; }
    void setConfig(const Config& config) { config_ = config; }
    
    // === Utility ===
    
    // Memory management
    void compact();
    double getFragmentation() const;
    
    // NOTE: Additional utility methods (time management, validation, serialization)
    // are not implemented in this stub version
    
private:
    // === Private Data Members ===
    
    // Templates
    std::map<int, FragmentTemplate> templates_;
    std::map<std::string, int> templateNameMap_;
    
    // Instances  
    std::map<int, FragmentInstance> instances_;
    std::set<int> activeInstances_;                     // Active instance IDs
    std::set<int> ghostInstances_;                      // Ghost instance IDs
    std::map<int, std::set<int>> templateInstances_;    // Instances by template
    std::map<int, std::queue<int>> ghostQueues_;        // Ghost queues by template
    
    // Active management
    std::shared_ptr<ActivePool> pool_;                  // Underlying storage
    
    // Configuration and statistics
    Config config_;
    mutable Statistics stats_;
    
    // Current state
    int nextInstanceId_ = 0;
    int nextTemplateId_ = 0;
    double currentStep_ = 0.0;
    
    // === Private Helper Functions ===
    
    // Instance management helpers
    int allocateInstanceSlot();
    void freeInstanceSlot(int slot);
    std::vector<MCAtom> transformAtoms(const std::vector<MCAtom>& atoms,
                                       const Vector3& position,
                                       const Quaternion& orientation) const;
    
    // Ghost management helpers  
    void convertToGhost(int instanceId);
    void recycleGhostSlot(int instanceId);
    int selectGhostForRecycling(int templateId);
    
    // Statistics helpers
    void updateStatistics(int templateId, const std::string& operation);
    void updateAverageLifetime(int templateId, double lifetime);
    void updateAcceptanceRate(int templateId, double rate);
    
    // Memory management helpers
    bool shouldCompact() const;
    void performCompaction();
    
    // Validation helpers
    bool isValidInstanceId(int id) const;
    bool isValidTemplateId(int id) const;
};

// ============================================================================
// Implementation of inline functions
// ============================================================================

inline bool FragmentReservoir::isValidInstanceId(int id) const {
    auto it = instances_.find(id);
    return it != instances_.end() && it->second.instanceId >= 0;
}

inline bool FragmentReservoir::isValidTemplateId(int id) const {
    return templates_.find(id) != templates_.end();
}

inline double FragmentReservoir::getTotalEnergy(int instanceId) const {
    auto it = instances_.find(instanceId);
    if (it == instances_.end()) return 0.0;
    return it->second.energy_total;
}

inline double FragmentReservoir::getVdwEnergy(int instanceId) const {
    auto it = instances_.find(instanceId);
    if (it == instances_.end()) return 0.0;
    return it->second.energy_vdw;
}

inline double FragmentReservoir::getElecEnergy(int instanceId) const {
    auto it = instances_.find(instanceId);
    if (it == instances_.end()) return 0.0;
    return it->second.energy_elec;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc