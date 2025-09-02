#include "CavityBiasCore.hpp"
#include "../common/MovementUtils.hpp"
#include <queue>
#include <algorithm>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

CavityBiasCore::CavityBiasCore(double gridSpacing, double probeRadius)
    : gridSpacing_(gridSpacing), probeRadius_(probeRadius), cacheValid_(false) {}

double CavityBiasCore::calculateCavityVolume(const MCState& state, CavityMode mode) {
    switch (mode) {
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
    if (!cacheValid_) {
        buildGrid(state);
    }
    
    switch (mode) {
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
    // Convert box from Angstroms to nm
    const double ANG_TO_NM = 0.1;
    grid_.box = Vector3(state.info.box[0] * ANG_TO_NM,
                       state.info.box[1] * ANG_TO_NM,
                       state.info.box[2] * ANG_TO_NM);
    
    // Setup grid dimensions
    grid_.nx = std::max(3, static_cast<int>(grid_.box.x / gridSpacing_));
    grid_.ny = std::max(3, static_cast<int>(grid_.box.y / gridSpacing_));
    grid_.nz = std::max(3, static_cast<int>(grid_.box.z / gridSpacing_));
    
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
    
    // Find clusters for CLUSTER_VOLUME mode
    clusters_ = findClusters();
    
    cacheValid_ = true;
}

void CavityBiasCore::markOccupied(const MCState& state) {
    const double ANG_TO_NM = 0.1;
    
    for (int resIdx = 0; resIdx < state.activeResidueCount; ++resIdx) {
        const MCResidue& res = state.residues[resIdx];
        if (!res.active) continue;
        
        for (int i = 0; i < res.atomCount; ++i) {
            int atomIdx = res.atomStart + i;
            if (atomIdx >= state.activeAtomCount) continue;
            
            const MCAtom& atom = state.atoms[atomIdx];
            Vector3 pos(atom.x * ANG_TO_NM, atom.y * ANG_TO_NM, atom.z * ANG_TO_NM);
            
            // Get effective radius (use LJ sigma if available)
            double radius = probeRadius_;
            if (atom.type < state.forcefield.numTotalTypes) {
                radius = 0.5 * state.forcefield.ljSigma[atom.type] * ANG_TO_NM + probeRadius_;
            }
            
            // Mark grid points within radius as occupied
            int iMin = std::max(0, (int)((pos.x - radius) / grid_.spacing.x));
            int iMax = std::min(grid_.nx - 1, (int)((pos.x + radius) / grid_.spacing.x));
            int jMin = std::max(0, (int)((pos.y - radius) / grid_.spacing.y));
            int jMax = std::min(grid_.ny - 1, (int)((pos.y + radius) / grid_.spacing.y));
            int kMin = std::max(0, (int)((pos.z - radius) / grid_.spacing.z));
            int kMax = std::min(grid_.nz - 1, (int)((pos.z + radius) / grid_.spacing.z));
            
            for (int i = iMin; i <= iMax; ++i) {
                for (int j = jMin; j <= jMax; ++j) {
                    for (int k = kMin; k <= kMax; ++k) {
                        Vector3 gridPos(i * grid_.spacing.x,
                                       j * grid_.spacing.y,
                                       k * grid_.spacing.z);
                        if (distance(pos, gridPos) < radius) {
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
    
    // Mode A: Simple cavity fraction
    int cavityCount = cavityPoints_.size();
    
    double voxelVolume = grid_.spacing.x * grid_.spacing.y * grid_.spacing.z;
    return cavityCount * voxelVolume;  // nm³
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
    if (cavityPoints_.empty()) {
        return Vector3(grid_.box.x * 0.5, grid_.box.y * 0.5, grid_.box.z * 0.5);
    }
    
    int idx = utils::RandomUtils::uniformInt(0, cavityPoints_.size() - 1);
    return cavityPoints_[idx];
}

Vector3 CavityBiasCore::sampleClusterVolume() {
    if (clusters_.empty()) {
        return sampleFastApprox();
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
    
    // Select random point within cluster
    const auto& cluster = clusters_[clusterIdx];
    int pointIdx = utils::RandomUtils::uniformInt(0, cluster.size() - 1);
    return cavityPoints_[cluster[pointIdx]];
}

Vector3 CavityBiasCore::sampleLocalVeff(const MCState& /* state */) {
    // For simplicity, use cluster sampling as base
    Vector3 basePos = sampleClusterVolume();
    
    // Add small random offset for continuous sampling
    double offset = 0.01;  // nm
    basePos.x += (utils::RandomUtils::uniform() - 0.5) * offset;
    basePos.y += (utils::RandomUtils::uniform() - 0.5) * offset;
    basePos.z += (utils::RandomUtils::uniform() - 0.5) * offset;
    
    // Apply PBC
    basePos.x = fmod(basePos.x + grid_.box.x, grid_.box.x);
    basePos.y = fmod(basePos.y + grid_.box.y, grid_.box.y);
    basePos.z = fmod(basePos.z + grid_.box.z, grid_.box.z);
    
    return basePos;
}

std::vector<std::vector<int>> CavityBiasCore::findClusters() {
    std::vector<std::vector<int>> clusters;
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
            
            // Check neighbors
            for (size_t j = 0; j < cavityPoints_.size(); ++j) {
                if (!visited[j] && 
                    distance(cavityPoints_[current], cavityPoints_[j]) < grid_.spacing.x * 1.5) {
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

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc