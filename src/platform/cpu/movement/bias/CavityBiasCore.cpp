#include "CavityBiasCore.hpp"
#include "../common/MovementUtils.hpp"
#include <queue>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

CavityBiasCore::CavityBiasCore(double gridSpacing, double probeRadius)
    : gridSpacing_(gridSpacing),
      probeRadius_(probeRadius),
      activeGridSpacing_(gridSpacing),
      activeProbeRadius_(probeRadius),
      activeSpeciesId_(-1),
      cacheValid_(false) {}

double CavityBiasCore::calculateCavityVolume(const MCState& state, CavityMode mode, int speciesId) {
    applySpeciesParameters(speciesId);
    // Auto-detect ideal gas and use appropriate mode
    CavityMode effectiveMode = mode;
    if (mode == CavityMode::CLUSTER_VOLUME || mode == CavityMode::LOCAL_VEFF) {
        if (isIdealGas(state)) {
            // For ideal gas, use FAST_APPROX to reduce discretization errors
            effectiveMode = CavityMode::FAST_APPROX;
            // Debug: log mode switch (disabled)
            // std::cerr << "CavityBiasCore: Ideal gas detected, switching from mode "
            //           << static_cast<int>(mode) << " to FAST_APPROX" << std::endl;
        }
    }

    switch (effectiveMode) {
        case CavityMode::FAST_APPROX:
            return calculateFastApprox(state);
        case CavityMode::CLUSTER_VOLUME:
            return calculateClusterVolume(state);
        case CavityMode::LOCAL_VEFF:
            // For global volume, use cluster volume as approximation
            return calculateClusterVolume(state);
        default:
            return calculateFastApprox(state);
    }
}

Vector3 CavityBiasCore::proposeCavityPosition(const MCState& state, CavityMode mode) {
    applySpeciesParameters(activeSpeciesId_);
    if (!cacheValid_) {
        buildGrid(state);
    }

    // Auto-detect ideal gas and use appropriate mode
    CavityMode effectiveMode = mode;
    if (mode == CavityMode::CLUSTER_VOLUME || mode == CavityMode::LOCAL_VEFF) {
        if (isIdealGas(state)) {
            effectiveMode = CavityMode::FAST_APPROX;
            // std::cerr << "CavityBiasCore::proposeCavityPosition: Ideal gas, using FAST_APPROX" << std::endl;
        }
    }

    switch (effectiveMode) {
        case CavityMode::FAST_APPROX:
            return sampleFastApprox();
        case CavityMode::CLUSTER_VOLUME:
            return sampleClusterVolume();
        case CavityMode::LOCAL_VEFF:
            return sampleLocalVeff(state);
        default:
            return sampleFastApprox();
    }
}

void CavityBiasCore::buildGrid(const MCState& state) {
    // Validate box dimensions
    if (state.info.box[0] <= 0 || state.info.box[1] <= 0 || state.info.box[2] <= 0) {
        // Invalid box dimensions - clear cache and return
        cacheValid_ = false;
        cavityPoints_.clear();
        clusters_.clear();
        grid_.occupied.clear();
        return;
    }

    // Ensure grid spacing is positive
    double spacing = activeGridSpacing_ > 0 ? activeGridSpacing_ : gridSpacing_;
    if (spacing <= 0.0) {
        spacing = 0.25;  // Default to 0.25 nm
    }

    // Box is already in nm in MCState
    grid_.box = Vector3(state.info.box[0],
                       state.info.box[1],
                       state.info.box[2]);

    // Setup grid dimensions - use round for stability
    grid_.nx = std::max(3, static_cast<int>(std::round(grid_.box.x / spacing)));
    grid_.ny = std::max(3, static_cast<int>(std::round(grid_.box.y / spacing)));
    grid_.nz = std::max(3, static_cast<int>(std::round(grid_.box.z / spacing)));

    grid_.origin = Vector3(0, 0, 0);
    grid_.spacing = Vector3(grid_.box.x / grid_.nx,
                           grid_.box.y / grid_.ny,
                           grid_.box.z / grid_.nz);

    // Initialize occupancy
    int totalPoints = grid_.nx * grid_.ny * grid_.nz;
    grid_.occupied.assign(totalPoints, false);

    // Mark occupied regions
    markOccupied(state);

    // Build cavity point list
    cavityPoints_.clear();
    clusters_.clear();

    // For performance, limit cavity point collection for very large grids
    bool buildFullList = (totalPoints <= 50000);

    if (buildFullList) {
        // Build complete cavity point list for smaller grids
        for (int i = 0; i < grid_.nx; ++i) {
            for (int j = 0; j < grid_.ny; ++j) {
                for (int k = 0; k < grid_.nz; ++k) {
                    int idx = grid_.getIndex(i, j, k);
                    if (!grid_.occupied[idx]) {
                        Vector3 pos(grid_.origin.x + i * grid_.spacing.x,
                                   grid_.origin.y + j * grid_.spacing.y,
                                   grid_.origin.z + k * grid_.spacing.z);
                        cavityPoints_.push_back(pos);
                    }
                }
            }
        }

        // Find clusters only if cavity points list is reasonable size
        if (cavityPoints_.size() < 5000) {
            clusters_ = findClusters();
        } else {
            // For large cavity sets, treat as single cluster
            std::vector<int> singleCluster;
            for (size_t i = 0; i < cavityPoints_.size(); ++i) {
                singleCluster.push_back(i);
            }
            clusters_.push_back(singleCluster);
        }
    } else {
        // For very large grids, sample a subset of cavity points
        int sampleStride = std::max(2, static_cast<int>(std::cbrt(totalPoints / 10000.0)));

        for (int i = 0; i < grid_.nx; i += sampleStride) {
            for (int j = 0; j < grid_.ny; j += sampleStride) {
                for (int k = 0; k < grid_.nz; k += sampleStride) {
                    int idx = grid_.getIndex(i, j, k);
                    if (!grid_.occupied[idx]) {
                        Vector3 pos(grid_.origin.x + i * grid_.spacing.x,
                                   grid_.origin.y + j * grid_.spacing.y,
                                   grid_.origin.z + k * grid_.spacing.z);
                        cavityPoints_.push_back(pos);
                    }
                }
            }
        }

        // Treat sampled points as single cluster
        if (!cavityPoints_.empty()) {
            std::vector<int> singleCluster;
            for (size_t i = 0; i < cavityPoints_.size(); ++i) {
                singleCluster.push_back(i);
            }
            clusters_.push_back(singleCluster);
        }
    }

    cacheValid_ = true;
}

void CavityBiasCore::markOccupied(const MCState& state) {
    // Coordinates and sigma are already in nm

    // Ensure we don't exceed residues vector size
    int maxResIdx = std::min(state.activeResidueCount, static_cast<int>(state.residues.size()));

    for (int resIdx = 0; resIdx < maxResIdx; ++resIdx) {
        const MCResidue& res = state.residues[resIdx];
        if (!res.active) continue;

        for (int i = 0; i < res.atomCount; ++i) {
            int atomIdx = res.atomStart + i;
            // Check for negative index (e.g., when atomStart is -1)
            // Also ensure we don't exceed atoms vector size
            if (atomIdx < 0 || atomIdx >= state.activeAtomCount ||
                atomIdx >= static_cast<int>(state.atoms.size())) {
                continue;
            }

            const MCAtom& atom = state.atoms[atomIdx];
            Vector3 pos(atom.x, atom.y, atom.z);  // Already in nm

            // Get effective radius (use LJ sigma if available)
            double radius = activeProbeRadius_;
            if (atom.type < state.forcefield.numTotalTypes &&
                !state.forcefield.ljSigma.empty()) {
                // ljSigma is an NxN matrix, get diagonal element
                int idx = atom.type * state.forcefield.numTotalTypes + atom.type;
                if (idx < static_cast<int>(state.forcefield.ljSigma.size())) {
                    radius = 0.5 * state.forcefield.ljSigma[idx] + probeRadius_;  // Already in nm
                }
            }

            // Mark grid points within radius as occupied (with PBC)
            int iMin = (int)((pos.x - radius) / grid_.spacing.x) - 1;
            int iMax = (int)((pos.x + radius) / grid_.spacing.x) + 1;
            int jMin = (int)((pos.y - radius) / grid_.spacing.y) - 1;
            int jMax = (int)((pos.y + radius) / grid_.spacing.y) + 1;
            int kMin = (int)((pos.z - radius) / grid_.spacing.z) - 1;
            int kMax = (int)((pos.z + radius) / grid_.spacing.z) + 1;

            for (int ii = iMin; ii <= iMax; ++ii) {
                for (int jj = jMin; jj <= jMax; ++jj) {
                    for (int kk = kMin; kk <= kMax; ++kk) {
                        // Apply PBC wrapping
                        int i = ((ii % grid_.nx) + grid_.nx) % grid_.nx;
                        int j = ((jj % grid_.ny) + grid_.ny) % grid_.ny;
                        int k = ((kk % grid_.nz) + grid_.nz) % grid_.nz;

                        Vector3 gridPos(i * grid_.spacing.x,
                                       j * grid_.spacing.y,
                                       k * grid_.spacing.z);

                        // Calculate minimum image distance
                        double dx = gridPos.x - pos.x;
                        double dy = gridPos.y - pos.y;
                        double dz = gridPos.z - pos.z;

                        dx -= std::round(dx / grid_.box.x) * grid_.box.x;
                        dy -= std::round(dy / grid_.box.y) * grid_.box.y;
                        dz -= std::round(dz / grid_.box.z) * grid_.box.z;

                        double dist2 = dx*dx + dy*dy + dz*dz;
                        if (dist2 < radius * radius) {
                            grid_.occupied[grid_.getIndex(i, j, k)] = true;
                        }
                    }
                }
            }
        }
    }
}

double CavityBiasCore::calculateFastApprox(const MCState& state) {
    if (!cacheValid_) {
        buildGrid(state);
    }

    // Mode A: Simple cavity fraction (skip expensive cluster analysis)
    // For FAST mode, we don't need cavityPoints_ list, just count unoccupied voxels
    int cavityCount = 0;
    int totalCount = grid_.occupied.size();
    for (size_t i = 0; i < grid_.occupied.size(); ++i) {
        if (!grid_.occupied[i]) cavityCount++;
    }

    double voxelVolume = grid_.spacing.x * grid_.spacing.y * grid_.spacing.z;
    double cavityVolume = cavityCount * voxelVolume;  // nm³

    // Debug output (remove in production)
    if (std::getenv("DEBUG_CAVITY")) {
        std::cout << "[CavityBiasCore] FastApprox: cavity=" << cavityCount
                  << "/" << totalCount << " voxels, volume=" << cavityVolume
                  << " nm³, n_particles=" << state.activeResidueCount << std::endl;
    }

    return cavityVolume;
}

double CavityBiasCore::calculateClusterVolume(const MCState& state) {
    if (!cacheValid_) {
        buildGrid(state);
    }

    // Mode B: Sum of cluster volumes
    double totalVolume = 0.0;
    double voxelVolume = grid_.spacing.x * grid_.spacing.y * grid_.spacing.z;

    for (const auto& cluster : clusters_) {
        totalVolume += cluster.size() * voxelVolume;
    }

    return totalVolume;  // nm³
}

double CavityBiasCore::calculateLocalVeff(const MCState& /* state */, const Vector3& pos) {
    // Mode C: Local effective volume estimation
    const int subsamples = 100;
    const double subsampleRadius = 0.1;  // nm
    int accessibleCount = 0;

    for (int i = 0; i < subsamples; ++i) {
        // Random offset within subsample sphere
        double r = utils::RandomUtils::uniform() * subsampleRadius;
        double theta = utils::RandomUtils::uniform() * 2 * M_PI;
        double phi = std::acos(2 * utils::RandomUtils::uniform() - 1);

        Vector3 offset(r * sin(phi) * cos(theta),
                      r * sin(phi) * sin(theta),
                      r * cos(phi));

        Vector3 samplePos = pos + offset;

        // Check if sample position is accessible
        int gi = static_cast<int>(samplePos.x / grid_.spacing.x);
        int gj = static_cast<int>(samplePos.y / grid_.spacing.y);
        int gk = static_cast<int>(samplePos.z / grid_.spacing.z);

        if (grid_.isValid(gi, gj, gk) && !grid_.occupied[grid_.getIndex(gi, gj, gk)]) {
            accessibleCount++;
        }
    }

    double sphereVolume = (4.0/3.0) * M_PI * subsampleRadius * subsampleRadius * subsampleRadius;
    return (accessibleCount / static_cast<double>(subsamples)) * sphereVolume;
}

Vector3 CavityBiasCore::sampleFastApprox() {
    // If we have cavity points, use them
    if (!cavityPoints_.empty()) {
        int idx = utils::RandomUtils::uniformInt(0, cavityPoints_.size() - 1);
        return cavityPoints_[idx];
    }

    // Otherwise, directly sample from unoccupied grid points
    // This is slower but ensures we always return a valid position
    std::vector<int> unoccupiedIndices;
    for (int i = 0; i < grid_.nx; ++i) {
        for (int j = 0; j < grid_.ny; ++j) {
            for (int k = 0; k < grid_.nz; ++k) {
                int idx = grid_.getIndex(i, j, k);
                if (!grid_.occupied[idx]) {
                    unoccupiedIndices.push_back(idx);
                }
            }
        }
    }

    if (unoccupiedIndices.empty()) {
        // No cavity available, return center
        return Vector3(grid_.box.x * 0.5, grid_.box.y * 0.5, grid_.box.z * 0.5);
    }

    // Randomly select an unoccupied voxel
    int selectedIdx = unoccupiedIndices[utils::RandomUtils::uniformInt(0, unoccupiedIndices.size() - 1)];

    // Convert grid index back to i,j,k
    int k = selectedIdx / (grid_.nx * grid_.ny);
    int j = (selectedIdx % (grid_.nx * grid_.ny)) / grid_.nx;
    int i = selectedIdx % grid_.nx;

    // Return center of the voxel with small random offset for better sampling
    return Vector3(
        grid_.origin.x + (i + 0.5) * grid_.spacing.x,
        grid_.origin.y + (j + 0.5) * grid_.spacing.y,
        grid_.origin.z + (k + 0.5) * grid_.spacing.z
    );
}

Vector3 CavityBiasCore::sampleClusterVolume() {
    // Ensure we have valid cavity points
    if (cavityPoints_.empty()) {
        // Fallback to direct grid sampling
        return sampleFastApprox();
    }

    if (clusters_.empty()) {
        // No clusters, but have cavity points - sample uniformly
        int idx = utils::RandomUtils::uniformInt(0, cavityPoints_.size() - 1);
        return cavityPoints_[idx];
    }

    // Build cumulative distribution based on cluster volumes
    std::vector<double> cumulative(clusters_.size() + 1, 0.0);
    for (size_t i = 0; i < clusters_.size(); ++i) {
        cumulative[i + 1] = cumulative[i] + clusters_[i].size();
    }

    // Select cluster weighted by volume
    double r = utils::RandomUtils::uniform() * cumulative.back();
    auto it = std::upper_bound(cumulative.begin(), cumulative.end(), r);
    int clusterIdx = std::distance(cumulative.begin(), it) - 1;

    // Bounds check
    if (clusterIdx < 0 || clusterIdx >= static_cast<int>(clusters_.size())) {
        return sampleFastApprox();
    }

    // Select random point within cluster
    const auto& cluster = clusters_[clusterIdx];
    if (cluster.empty()) {
        return sampleFastApprox();
    }

    int pointIdx = utils::RandomUtils::uniformInt(0, cluster.size() - 1);

    // Bounds check for cavity point access
    if (cluster[pointIdx] >= 0 && cluster[pointIdx] < static_cast<int>(cavityPoints_.size())) {
        return cavityPoints_[cluster[pointIdx]];
    }

    // Fallback
    return sampleFastApprox();
}

Vector3 CavityBiasCore::sampleLocalVeff(const MCState& /* state */) {
    // Start with cluster sampling
    Vector3 basePos = sampleClusterVolume();

    // Try to add small random offset while staying in cavity
    const double offset = 0.01;  // nm
    const int maxAttempts = 10;

    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        Vector3 candidate = basePos;
        candidate.x += (utils::RandomUtils::uniform() - 0.5) * offset;
        candidate.y += (utils::RandomUtils::uniform() - 0.5) * offset;
        candidate.z += (utils::RandomUtils::uniform() - 0.5) * offset;

        // Apply PBC
        while (candidate.x < 0) candidate.x += grid_.box.x;
        while (candidate.x >= grid_.box.x) candidate.x -= grid_.box.x;
        while (candidate.y < 0) candidate.y += grid_.box.y;
        while (candidate.y >= grid_.box.y) candidate.y -= grid_.box.y;
        while (candidate.z < 0) candidate.z += grid_.box.z;
        while (candidate.z >= grid_.box.z) candidate.z -= grid_.box.z;

        // Check if still in cavity
        int gi = static_cast<int>(candidate.x / grid_.spacing.x);
        int gj = static_cast<int>(candidate.y / grid_.spacing.y);
        int gk = static_cast<int>(candidate.z / grid_.spacing.z);

        if (grid_.isValid(gi, gj, gk) && !grid_.occupied[grid_.getIndex(gi, gj, gk)]) {
            return candidate;  // Found valid position in cavity
        }
    }

    // If no valid offset found, return original position
    return basePos;
}

std::vector<std::vector<int>> CavityBiasCore::findClusters() {
    std::vector<std::vector<int>> clusters;

    // Skip clustering for large cavity sets (performance optimization)
    if (cavityPoints_.size() > 2000) {
        // Treat entire cavity as single cluster
        std::vector<int> singleCluster;
        for (size_t i = 0; i < cavityPoints_.size(); ++i) {
            singleCluster.push_back(i);
        }
        clusters.push_back(singleCluster);
        return clusters;
    }

    std::vector<bool> visited(cavityPoints_.size(), false);

    for (size_t i = 0; i < cavityPoints_.size(); ++i) {
        if (visited[i]) continue;

        std::vector<int> cluster;
        std::queue<int> queue;
        queue.push(i);
        visited[i] = true;

        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();
            cluster.push_back(current);

            // Optimization: only check nearby points using grid locality
            // Instead of O(n²), limit search to reasonable neighbors
            const Vector3& currentPos = cavityPoints_[current];

            for (size_t j = 0; j < cavityPoints_.size(); ++j) {
                if (visited[j]) continue;

                const Vector3& candidatePos = cavityPoints_[j];

                // Quick rejection based on Manhattan distance
                double dx = std::abs(candidatePos.x - currentPos.x);
                double dy = std::abs(candidatePos.y - currentPos.y);
                double dz = std::abs(candidatePos.z - currentPos.z);

                // Apply PBC for Manhattan distance
                dx = std::min(dx, grid_.box.x - dx);
                dy = std::min(dy, grid_.box.y - dy);
                dz = std::min(dz, grid_.box.z - dz);

                // Skip if Manhattan distance is too large
                if (dx > grid_.spacing.x * 2 || dy > grid_.spacing.y * 2 || dz > grid_.spacing.z * 2) {
                    continue;
                }

                // Detailed distance check
                if (distance(currentPos, candidatePos) < grid_.spacing.x * 1.5) {
                    queue.push(j);
                    visited[j] = true;
                }
            }
        }

        clusters.push_back(cluster);
    }

    return clusters;
}

double CavityBiasCore::distance(const Vector3& a, const Vector3& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;

    // Apply PBC
    dx = dx - grid_.box.x * round(dx / grid_.box.x);
    dy = dy - grid_.box.y * round(dy / grid_.box.y);
    dz = dz - grid_.box.z * round(dz / grid_.box.z);

    return sqrt(dx*dx + dy*dy + dz*dz);
}

bool CavityBiasCore::isIdealGas(const MCState& state) const {
    // Check if all LJ epsilon values are essentially zero
    // This indicates an ideal gas system with no intermolecular interactions
    const double eps_threshold = 1e-6;  // Threshold for "zero" epsilon

    const int n = state.forcefield.numTotalTypes;
    if (n <= 0) return true;  // No types defined, treat as ideal gas

    const auto& eps = state.forcefield.ljEps;
    // Check if ljEps is properly sized for NxN matrix
    if (eps.size() < static_cast<size_t>(n * n)) {
        // ljEps not properly initialized, treat as ideal gas for safety
        return true;
    }

    // Check diagonal elements (self-interactions) of the epsilon matrix
    bool isIdeal = true;
    for (int i = 0; i < n; ++i) {
        // ljEps is an NxN matrix, get diagonal element
        int idx = i * n + i;
        if (eps[idx] > eps_threshold) {
            isIdeal = false;
            break;
        }
    }

    // Debug output (disabled)
    // std::cerr << "CavityBiasCore::isIdealGas: numTypes=" << n
    //           << ", ljEps.size=" << eps.size()
    //           << ", result=" << (isIdeal ? "true" : "false") << std::endl;

    return isIdeal;
}

void CavityBiasCore::setSpeciesParameters(int speciesId, double gridSpacing, double probeRadius, int maskId) {
    if (speciesId < 0) {
        return;
    }
    SpeciesParameters cfg;
    cfg.gridSpacingNm = gridSpacing;
    cfg.probeRadiusNm = probeRadius;
    cfg.maskId = maskId;
    speciesParams_[speciesId] = cfg;
    cacheValid_ = false;
}

void CavityBiasCore::clearSpeciesParameters() {
    speciesParams_.clear();
    activeSpeciesId_ = -1;
    activeGridSpacing_ = gridSpacing_;
    activeProbeRadius_ = probeRadius_;
    cacheValid_ = false;
}

CavityBiasCore::SpeciesParameters CavityBiasCore::resolveSpeciesParameters(int speciesId) const {
    auto it = speciesParams_.find(speciesId);
    if (it != speciesParams_.end()) {
        return it->second;
    }
    return {};
}

void CavityBiasCore::applySpeciesParameters(int speciesId) {
    double spacing = gridSpacing_;
    double probe = probeRadius_;
    if (speciesId >= 0) {
        SpeciesParameters cfg = resolveSpeciesParameters(speciesId);
        if (cfg.gridSpacingNm > 0.0) {
            spacing = cfg.gridSpacingNm;
        }
        if (cfg.probeRadiusNm > 0.0) {
            probe = cfg.probeRadiusNm;
        }
    }

    const bool changed = (speciesId != activeSpeciesId_) ||
                         (std::abs(spacing - activeGridSpacing_) > 1e-12) ||
                         (std::abs(probe - activeProbeRadius_) > 1e-12);
    if (changed) {
        activeSpeciesId_ = speciesId;
        activeGridSpacing_ = spacing;
        activeProbeRadius_ = probe;
        cacheValid_ = false;
    }
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
