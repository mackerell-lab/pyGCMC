#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_CAVITY_BIAS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_CAVITY_BIAS_HPP

#include <vector>
#include <memory>
#include <set>
#include <mutex>
#include "../common/MovementUtils.hpp"

// Forward declarations
namespace pygcmc {
namespace model {
namespace montecarlo {
    struct MCState;
    struct MCAtom;
}
}
}

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCState;
using model::montecarlo::MCAtom;

/**
 * Cavity Manager for identifying and managing cavities in the system
 * Used for cavity-biased insertion to improve acceptance rates
 */
class CavityManager {
public:
    // Grid structure for cavity identification
    struct Grid3D {
        std::vector<bool> occupied;      // Occupancy grid
        Vector3 origin;                  // Grid origin
        Vector3 spacing;                 // Grid spacing
        Vector3 boxSize;                 // System box size
        int nx, ny, nz;                  // Grid dimensions
        
        // Get linear index from grid coordinates
        int getIndex(int i, int j, int k) const {
            return i + nx * (j + ny * k);
        }
        
        // Get grid coordinates from linear index
        void getCoords(int index, int& i, int& j, int& k) const {
            k = index / (nx * ny);
            j = (index % (nx * ny)) / nx;
            i = index % nx;
        }
        
        // Check if grid point is valid
        bool isValid(int i, int j, int k) const {
            return i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz;
        }
    };
    
    // Constructor
    explicit CavityManager(double gridSpacing = 2.0, double probeRadius = 1.4);
    
    // Destructor
    ~CavityManager();
    
    // Main cavity finding function
    std::vector<Vector3> findCavities(const MCState& state);
    
    // Calculate cavity bias factor for acceptance probability
    double calculateCavityBiasFactor(const MCState& state);
    
    // Check if a position is in a cavity
    bool isInCavity(const Vector3& position, const MCState& state);
    
    // Cache management
    void invalidateCache();
    bool isCacheValid() const { return cacheValid_; }
    
    // Configuration
    void setGridSpacing(double spacing) { gridSpacing_ = spacing; invalidateCache(); }
    void setProbeRadius(double radius) { probeRadius_ = radius; invalidateCache(); }
    double getGridSpacing() const { return gridSpacing_; }
    double getProbeRadius() const { return probeRadius_; }
    
    // Grid access
    const Grid3D& getGrid() const { return grid_; }
    int getTotalGridPoints() const { return grid_.nx * grid_.ny * grid_.nz; }
    int getCavityCount() const { return static_cast<int>(cavityCache_.size()); }
    
    // Statistics (P2 enhanced)
    struct Statistics {
        // Basic grid statistics
        int totalGridPoints = 0;
        int occupiedPoints = 0;
        int cavityPoints = 0;
        double occupancyRatio = 0.0;
        double cavityRatio = 0.0;
        
        // Cache performance
        int cacheHits = 0;
        int cacheMisses = 0;
        
        // Cluster analysis
        int clusterCount = 0;
        int largestClusterSize = 0;
        int averageClusterSize = 0;
        
        // P2: Enhanced metrics
        double buildTimeMs = 0.0;           // Time to build cavity grid (ms)
        double lastFindTimeMs = 0.0;        // Last cavity finding time (ms)
        int colorClassCount = 0;            // Number of color classes
        
        // Cavity distribution per color class
        struct ColorClassStats {
            int minCavities = 0;
            int medianCavities = 0;
            int p95Cavities = 0;
            double avgCavities = 0.0;
        } colorClassStats;
        
        // Incremental update stats (P3 prep)
        int incrementalUpdates = 0;
        int fullRebuilds = 0;
        double dirtyRatio = 0.0;            // Fraction of grid marked dirty
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
    // Advanced cavity identification with clustering
    struct CavityCluster {
        std::vector<Vector3> positions;
        Vector3 center;
        double volume;
        int id;
    };
    
    std::vector<CavityCluster> findCavityClusters(const MCState& state);
    
    // New clustering methods
    std::vector<CavityCluster> findCavityClustersFloodFill(const MCState& state);
    Vector3 selectFromCluster(const CavityCluster& cluster);
    std::vector<Vector3> getClusterCenters(const MCState& state, int maxClusters = 100);
    
    // Color-based independent selection
    std::vector<Vector3> selectIndependentCavities(const MCState& state, 
                                                   double minSeparationNm,
                                                   int maxPoints = 100);
    
    // Automatic mode selection
    bool shouldUseClustering(const MCState& state);
    
private:
    // Configuration
    double gridSpacing_;                    // Grid spacing in Angstroms
    double probeRadius_;                    // Probe radius for cavity detection
    
    // Grid structure
    Grid3D grid_;
    
    // Cache
    std::vector<Vector3> cavityCache_;     // Cached cavity positions
    bool cacheValid_;                      // Whether cache is valid
    Vector3 lastBoxSize_;                  // Last box size for cache invalidation (per-instance)
    
    // Statistics
    mutable Statistics stats_;
    
    // Thread safety
    mutable std::mutex cacheMutex_;        // Protects cache and lastBoxSize
    
    // Helper functions
    void initializeGrid(const Vector3& boxSize);
    void markOccupiedRegions(const MCState& state);
    void markOccupiedRegion(const Vector3& center, double radius);
    Vector3 gridToPosition(int i, int j, int k) const;
    bool checkCavity(const Vector3& position, const MCState& state) const;
    
    // Distance calculation with PBC
    double distance(const Vector3& pos1, const Vector3& pos2) const;
    
    // Clustering helper
    void clusterCavities(std::vector<CavityCluster>& clusters);
    int findCluster(const Vector3& pos, const std::vector<CavityCluster>& clusters, double threshold);
};

/**
 * Cavity-biased insertion manager
 */
class CavityBiasInsertion {
public:
    // Constructor
    explicit CavityBiasInsertion(CavityManager* cavityManager);
    
    // Select insertion position with cavity bias
    Vector3 selectInsertionPosition(const MCState& state, bool& usedCavity);
    
    // Calculate acceptance probability with cavity bias
    double calculateAcceptanceProbability(
        int n,                      // Current number of molecules
        double deltaE,              // Energy change
        double beta,                // 1/kT
        double chemPotential,       // Chemical potential
        double volumeNm3,           // System volume
        bool usedCavity             // Whether cavity was used
    );
    
    // Configuration
    void setUseCavityBias(bool use) { useCavityBias_ = use; }
    bool getUseCavityBias() const { return useCavityBias_; }
    
    // Statistics
    struct Statistics {
        int totalInsertions = 0;
        int cavityInsertions = 0;
        int randomInsertions = 0;
        double averageCavityBias = 0.0;
        double averageAcceptance = 0.0;
    };
    
    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics();
    
private:
    CavityManager* cavityManager_;
    bool useCavityBias_;
    Statistics stats_;
    
    // Helper functions
    Vector3 selectRandomPosition(const Vector3& boxSize);
    void updateStatistics(bool usedCavity, double cavityBias, double acceptance);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_CAVITY_BIAS_HPP