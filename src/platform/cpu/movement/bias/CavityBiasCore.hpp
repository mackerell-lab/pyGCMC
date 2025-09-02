#ifndef CAVITY_BIAS_CORE_HPP
#define CAVITY_BIAS_CORE_HPP

#include <vector>
#include <memory>
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../common/MovementUtils.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCState;

// Cavity bias modes based on three-tier strategy
enum class CavityMode {
    FAST_APPROX = 0,    // Global cavity fraction (default)
    CLUSTER_VOLUME = 1,  // Cluster-based accurate sampling
    LOCAL_VEFF = 2       // Local effective volume (highest accuracy)
};

// Minimal cavity grid structure
struct CavityGrid {
    Vector3 origin;
    Vector3 spacing;  // Grid spacing in nm
    Vector3 box;      // Box dimensions in nm
    int nx, ny, nz;   // Grid dimensions
    std::vector<bool> occupied;
    
    int getIndex(int i, int j, int k) const {
        return i + j * nx + k * nx * ny;
    }
    
    bool isValid(int i, int j, int k) const {
        return i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz;
    }
};

// Core cavity bias calculator
class CavityBiasCore {
public:
    CavityBiasCore(double gridSpacing = 0.25, double probeRadius = 0.14);
    ~CavityBiasCore() = default;
    
    // Main interface
    double calculateCavityVolume(const MCState& state, CavityMode mode);
    Vector3 proposeCavityPosition(const MCState& state, CavityMode mode);
    void invalidateCache() { cacheValid_ = false; }
    
    // Configuration
    void setGridSpacing(double spacing) { gridSpacing_ = spacing; }
    void setProbeRadius(double radius) { probeRadius_ = radius; }
    
private:
    // Grid building
    void buildGrid(const MCState& state);
    void markOccupied(const MCState& state);
    
    // Mode-specific calculations
    double calculateFastApprox(const MCState& state);
    double calculateClusterVolume(const MCState& state);
    double calculateLocalVeff(const MCState& state, const Vector3& pos);
    
    // Sampling
    Vector3 sampleFastApprox();
    Vector3 sampleClusterVolume();
    Vector3 sampleLocalVeff(const MCState& state);
    
    // Helper functions
    std::vector<std::vector<int>> findClusters();
    double distance(const Vector3& a, const Vector3& b) const;
    
    // Member variables
    double gridSpacing_;  // nm
    double probeRadius_;  // nm
    CavityGrid grid_;
    std::vector<Vector3> cavityPoints_;
    std::vector<std::vector<int>> clusters_;
    bool cacheValid_;
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // CAVITY_BIAS_CORE_HPP