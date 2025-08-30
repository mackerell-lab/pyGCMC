#include "CavityBias.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../common/MovementUtils.hpp"
#include <cmath>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

// Constants for unit conversion
static constexpr double ANGSTROM_TO_NM = 0.1;
static constexpr double NM_TO_ANGSTROM = 10.0;

CavityManager::CavityManager(double gridSpacing, double probeRadius)
    : gridSpacing_(gridSpacing),      // in Angstroms
      probeRadius_(probeRadius),      // in Angstroms
      cacheValid_(false),
      lastBoxSize_(-1, -1, -1) {      // Initialize to invalid size
    resetStatistics();
}

CavityManager::~CavityManager() = default;

std::vector<Vector3> CavityManager::findCavities(const MCState& state) {
    // Note: Thread safety issue if called concurrently - would need mutex
    // std::lock_guard<std::mutex> guard(cacheMutex_);
    
    // Get box dimensions from MCState (already in nm!)
    Vector3 boxSizeNm(state.info.box[0], 
                      state.info.box[1], 
                      state.info.box[2]);
    
    // Check if cache is valid (also check if box size changed)
    bool boxChanged = (std::abs(boxSizeNm.x - lastBoxSize_.x) > 1e-6 ||
                      std::abs(boxSizeNm.y - lastBoxSize_.y) > 1e-6 ||
                      std::abs(boxSizeNm.z - lastBoxSize_.z) > 1e-6);
    lastBoxSize_ = boxSizeNm;
    
    if (cacheValid_ && !cavityCache_.empty() && !boxChanged) {
        stats_.cacheHits++;
        return cavityCache_;
    }
    
    stats_.cacheMisses++;
    
    // Initialize grid (expects nm)
    initializeGrid(boxSizeNm);
    
    // Mark occupied regions
    markOccupiedRegions(state);
    
    // Find cavity points
    cavityCache_.clear();
    cavityCache_.reserve(grid_.nx * grid_.ny * grid_.nz / 10);  // Estimate ~10% cavities
    
    for (int i = 0; i < grid_.nx; ++i) {
        for (int j = 0; j < grid_.ny; ++j) {
            for (int k = 0; k < grid_.nz; ++k) {
                int idx = grid_.getIndex(i, j, k);
                if (!grid_.occupied[idx]) {
                    // store cavity position in nm
                    Vector3 posNm = gridToPosition(i, j, k);
                    cavityCache_.push_back(posNm);
                }
            }
        }
    }
    
    // Update statistics
    stats_.totalGridPoints = grid_.nx * grid_.ny * grid_.nz;
    stats_.cavityPoints = static_cast<int>(cavityCache_.size());
    stats_.occupiedPoints = stats_.totalGridPoints - stats_.cavityPoints;
    stats_.cavityRatio = static_cast<double>(stats_.cavityPoints) / stats_.totalGridPoints;
    stats_.occupancyRatio = static_cast<double>(stats_.occupiedPoints) / stats_.totalGridPoints;
    
    cacheValid_ = true;
    return cavityCache_;
}

double CavityManager::calculateCavityBiasFactor(const MCState& state) {
    auto cavities = findCavities(state);
    
    if (stats_.totalGridPoints == 0) {
        return 1.0;  // No bias if no grid
    }
    
    // Cavity bias factor f_n = N_cavities / N_total
    double factor = static_cast<double>(cavities.size()) / stats_.totalGridPoints;
    
    // Avoid zero bias factor
    if (factor < 1e-10) {
        factor = 1e-10;
    }
    
    return factor;
}

bool CavityManager::isInCavity(const Vector3& position, const MCState& /*state*/) {
    // Position is already in nm, same as grid
    
    int i = static_cast<int>((position.x - grid_.origin.x) / grid_.spacing.x);
    int j = static_cast<int>((position.y - grid_.origin.y) / grid_.spacing.y);
    int k = static_cast<int>((position.z - grid_.origin.z) / grid_.spacing.z);
    
    if (!grid_.isValid(i, j, k)) {
        return false;
    }
    
    return !grid_.occupied[grid_.getIndex(i, j, k)];
}

void CavityManager::invalidateCache() {
    // Note: Thread safety issue if called concurrently - would need mutex  
    // std::lock_guard<std::mutex> guard(cacheMutex_);
    cacheValid_ = false;
    cavityCache_.clear();
}

void CavityManager::resetStatistics() {
    stats_ = Statistics();
}

std::vector<CavityManager::CavityCluster> CavityManager::findCavityClusters(const MCState& state) {
    auto cavities = findCavities(state);
    std::vector<CavityCluster> clusters;
    
    if (cavities.empty()) {
        return clusters;
    }
    
    // Simple clustering based on distance threshold
    double clusterThreshold = grid_.spacing.x * 1.5;  // Already in nm
    
    for (const auto& cavity : cavities) {
        int clusterIdx = findCluster(cavity, clusters, clusterThreshold);
        
        if (clusterIdx >= 0) {
            // Add to existing cluster
            clusters[clusterIdx].positions.push_back(cavity);
            // Update center
            Vector3& center = clusters[clusterIdx].center;
            int n = static_cast<int>(clusters[clusterIdx].positions.size());
            center = center * ((n - 1.0) / n) + cavity * (1.0 / n);
        } else {
            // Create new cluster
            CavityCluster newCluster;
            newCluster.id = static_cast<int>(clusters.size());
            newCluster.positions.push_back(cavity);
            newCluster.center = cavity;
            newCluster.volume = std::pow(grid_.spacing.x, 3);  // Initial volume in nm^3
            clusters.push_back(newCluster);
        }
    }
    
    // Calculate cluster volumes
    for (auto& cluster : clusters) {
        cluster.volume = cluster.positions.size() * std::pow(grid_.spacing.x, 3);
    }
    
    return clusters;
}

// Private helper functions

void CavityManager::initializeGrid(const Vector3& boxSize) {
    // Box size is in nm, grid spacing is in Angstroms, convert to nm
    double spacingNm = gridSpacing_ * ANGSTROM_TO_NM;
    
    // Calculate grid dimensions
    grid_.nx = std::max(1, static_cast<int>(std::ceil(boxSize.x / spacingNm)));
    grid_.ny = std::max(1, static_cast<int>(std::ceil(boxSize.y / spacingNm)));
    grid_.nz = std::max(1, static_cast<int>(std::ceil(boxSize.z / spacingNm)));
    
    
    // Set grid properties (store in nm for consistency with MCState)
    grid_.origin = Vector3(0.0, 0.0, 0.0);
    grid_.spacing = Vector3(spacingNm, spacingNm, spacingNm);
    grid_.boxSize = boxSize;  // Already in nm
    
    // Initialize occupancy grid
    int totalPoints = grid_.nx * grid_.ny * grid_.nz;
    grid_.occupied.clear();
    grid_.occupied.resize(totalPoints, false);
}

void CavityManager::markOccupiedRegions(const MCState& state) {
    // Process all active atoms
    for (int resIdx = 0; resIdx < state.activeResidueCount; ++resIdx) {
        const MCResidue& residue = state.residues[resIdx];
        if (!residue.active) continue;
        
        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                const MCAtom& atom = state.atoms[atomIdx];
                // Position is in nm (same as box dimensions)
                Vector3 posNm(atom.x, atom.y, atom.z);
                
                // Use sigma from force field if available
                double radiusNm = 0.15;  // Default fallback in nm (1.5 Angstroms)
                
                if (atom.type >= 0 && 
                    atom.type < static_cast<int>(state.forcefield.ljSigma.size())) {
                    // LJ sigma is in nm (check forcefield units)
                    // Use sigma/2 as atomic radius
                    double sigma = state.forcefield.ljSigma[atom.type];
                    if (sigma > 0) {
                        radiusNm = 0.5 * sigma;
                    }
                }
                
                // Add probe radius (convert from Angstroms to nm)
                double totalRadius = radiusNm + probeRadius_ * ANGSTROM_TO_NM;
                markOccupiedRegion(posNm, totalRadius);
            }
        }
    }
}

void CavityManager::markOccupiedRegion(const Vector3& center, double radius) {
    // Center and radius are in nm
    int iMin = static_cast<int>(std::floor((center.x - radius - grid_.origin.x) / grid_.spacing.x));
    int iMax = static_cast<int>(std::floor((center.x + radius - grid_.origin.x) / grid_.spacing.x));
    int jMin = static_cast<int>(std::floor((center.y - radius - grid_.origin.y) / grid_.spacing.y));
    int jMax = static_cast<int>(std::floor((center.y + radius - grid_.origin.y) / grid_.spacing.y));
    int kMin = static_cast<int>(std::floor((center.z - radius - grid_.origin.z) / grid_.spacing.z));
    int kMax = static_cast<int>(std::floor((center.z + radius - grid_.origin.z) / grid_.spacing.z));

    const int nx = grid_.nx, ny = grid_.ny, nz = grid_.nz;
    const Vector3 boxNm = grid_.boxSize; // Actually in nm, not Angstroms
    const double radiusSq = radius * radius;

    for (int ii = iMin; ii <= iMax; ++ii) {
        int i = (ii % nx + nx) % nx;
        for (int jj = jMin; jj <= jMax; ++jj) {
            int j = (jj % ny + ny) % ny;
            for (int kk = kMin; kk <= kMax; ++kk) {
                int k = (kk % nz + nz) % nz;

                Vector3 gridPos = gridToPosition(i, j, k); // nm

                // minimum-image displacement in nm
                double dx = gridPos.x - center.x;
                double dy = gridPos.y - center.y;
                double dz = gridPos.z - center.z;
                dx -= std::round(dx / boxNm.x) * boxNm.x;
                dy -= std::round(dy / boxNm.y) * boxNm.y;
                dz -= std::round(dz / boxNm.z) * boxNm.z;

                double distSq = dx*dx + dy*dy + dz*dz;

                if (distSq <= radiusSq) {
                    grid_.occupied[grid_.getIndex(i, j, k)] = true;
                }
            }
        }
    }
}

Vector3 CavityManager::gridToPosition(int i, int j, int k) const {
    // Return position in nm
    return Vector3(
        grid_.origin.x + i * grid_.spacing.x,
        grid_.origin.y + j * grid_.spacing.y,
        grid_.origin.z + k * grid_.spacing.z
    );
}

bool CavityManager::checkCavity(const Vector3& position, const MCState& state) const {
    // Position is in nm
    // Simple check: ensure minimum distance from all atoms
    double minDist = probeRadius_ * ANGSTROM_TO_NM;  // Convert probe radius to nm
    double minDistSq = minDist * minDist;
    
    for (int resIdx = 0; resIdx < state.activeResidueCount; ++resIdx) {
        const MCResidue& residue = state.residues[resIdx];
        if (!residue.active) continue;
        
        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                const MCAtom& atom = state.atoms[atomIdx];
                // Calculate distance squared with PBC (box in Angstroms)
                Vector3 atomPos(atom.x, atom.y, atom.z);
                double dist = utils::PBCUtils::minimumImageDistance(
                    position,
                    atomPos,
                    Vector3(state.info.box[0], state.info.box[1], state.info.box[2])
                );
                double distSq = dist * dist;  // Square the distance for comparison
                
                if (distSq < minDistSq) {
                    return false;
                }
            }
        }
    }
    
    return true;
}

double CavityManager::distance(const Vector3& pos1, const Vector3& pos2) const {
    Vector3 diff = pos2 - pos1;
    return diff.norm();
}

int CavityManager::findCluster(const Vector3& pos, const std::vector<CavityCluster>& clusters, double threshold) {
    for (size_t i = 0; i < clusters.size(); ++i) {
        if (distance(pos, clusters[i].center) < threshold) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// CavityBiasInsertion implementation

CavityBiasInsertion::CavityBiasInsertion(CavityManager* cavityManager)
    : cavityManager_(cavityManager),
      useCavityBias_(true) {
    resetStatistics();
}

Vector3 CavityBiasInsertion::selectInsertionPosition(const MCState& state, bool& usedCavity) {
    usedCavity = false;
    
    if (useCavityBias_ && cavityManager_) {
        auto cavities = cavityManager_->findCavities(state);
        
        if (!cavities.empty()) {
            // Select random cavity
            int idx = utils::RandomUtils::uniformInt(0, static_cast<int>(cavities.size()) - 1);
            usedCavity = true;
            stats_.cavityInsertions++;
            return cavities[idx];
        }
    }
    
    // Fall back to random position (box already in nm)
    stats_.randomInsertions++;
    return selectRandomPosition(Vector3(state.info.box[0], 
                                       state.info.box[1], 
                                       state.info.box[2]));
}

double CavityBiasInsertion::calculateAcceptanceProbability(
    int n,
    double deltaE,
    double beta,
    double chemPotential,
    double volumeNm3,
    bool usedCavity) {
    
    double cavityBias = 1.0;
    
    if (usedCavity && cavityManager_) {
        // Note: This should be called with the state BEFORE insertion
        // to get the correct cavity bias factor
        cavityBias = cavityManager_->getCavityCount() / 
                    static_cast<double>(cavityManager_->getTotalGridPoints());
        
        if (cavityBias < 1e-10) {
            cavityBias = 1e-10;
        }
    }
    
    // Use log-space calculator for numerical stability
    double prob = utils::LogSpaceCalculator::calculateInsertionProbability(
        n, deltaE, beta, chemPotential, cavityBias, volumeNm3, true
    );
    
    // Update statistics
    stats_.totalInsertions++;
    stats_.averageCavityBias = (stats_.averageCavityBias * (stats_.totalInsertions - 1) + cavityBias) / 
                               stats_.totalInsertions;
    stats_.averageAcceptance = (stats_.averageAcceptance * (stats_.totalInsertions - 1) + prob) / 
                              stats_.totalInsertions;
    
    return prob;
}

void CavityBiasInsertion::resetStatistics() {
    stats_ = Statistics();
}

Vector3 CavityBiasInsertion::selectRandomPosition(const Vector3& boxSize) {
    // Box size is in nm
    return Vector3(
        utils::RandomUtils::uniform(0.0, boxSize.x),
        utils::RandomUtils::uniform(0.0, boxSize.y),
        utils::RandomUtils::uniform(0.0, boxSize.z)
    );
}

void CavityBiasInsertion::updateStatistics(bool usedCavity, double cavityBias, double acceptance) {
    if (usedCavity) {
        stats_.cavityInsertions++;
    } else {
        stats_.randomInsertions++;
    }
    
    stats_.totalInsertions++;
    stats_.averageCavityBias = (stats_.averageCavityBias * (stats_.totalInsertions - 1) + cavityBias) / 
                               stats_.totalInsertions;
    stats_.averageAcceptance = (stats_.averageAcceptance * (stats_.totalInsertions - 1) + acceptance) / 
                              stats_.totalInsertions;
}

// Clustering Implementation

// Helper struct for flood-fill
struct GridPoint {
    int i, j, k;
    
    GridPoint(int ii, int jj, int kk) : i(ii), j(jj), k(kk) {}
    
    bool operator==(const GridPoint& other) const {
        return i == other.i && j == other.j && k == other.k;
    }
};

// Hash function for GridPoint
struct GridPointHash {
    std::size_t operator()(const GridPoint& p) const {
        return std::hash<int>()(p.i) ^ 
               (std::hash<int>()(p.j) << 1) ^ 
               (std::hash<int>()(p.k) << 2);
    }
};

std::vector<CavityManager::CavityCluster> CavityManager::findCavityClustersFloodFill(const MCState& state) {
    // First, build cavity grid if needed
    if (!cacheValid_) {
        findCavities(state);
    }
    
    std::vector<CavityCluster> clusters;
    
    // Create visited grid
    std::vector<bool> visited(grid_.nx * grid_.ny * grid_.nz, false);
    
    // Flood-fill to find connected components
    for (int i = 0; i < grid_.nx; ++i) {
        for (int j = 0; j < grid_.ny; ++j) {
            for (int k = 0; k < grid_.nz; ++k) {
                int idx = grid_.getIndex(i, j, k);
                
                // Skip if occupied or already visited
                if (grid_.occupied[idx] || visited[idx]) {
                    continue;
                }
                
                // Start new cluster
                CavityCluster cluster;
                cluster.id = static_cast<int>(clusters.size());
                
                // Flood-fill using BFS
                std::queue<GridPoint> queue;
                queue.push(GridPoint(i, j, k));
                visited[idx] = true;
                
                Vector3 centerSum(0, 0, 0);
                int count = 0;
                
                while (!queue.empty()) {
                    GridPoint current = queue.front();
                    queue.pop();
                    
                    // Add to cluster
                    Vector3 posNm = gridToPosition(current.i, current.j, current.k);  // Already in nm
                    cluster.positions.push_back(posNm);
                    centerSum = centerSum + posNm;
                    count++;
                    
                    // Check 6-connected neighbors
                    const int di[] = {1, -1, 0, 0, 0, 0};
                    const int dj[] = {0, 0, 1, -1, 0, 0};
                    const int dk[] = {0, 0, 0, 0, 1, -1};
                    
                    for (int n = 0; n < 6; ++n) {
                        int ni = current.i + di[n];
                        int nj = current.j + dj[n];
                        int nk = current.k + dk[n];
                        
                        // Apply PBC
                        ni = (ni + grid_.nx) % grid_.nx;
                        nj = (nj + grid_.ny) % grid_.ny;
                        nk = (nk + grid_.nz) % grid_.nz;
                        
                        int nidx = grid_.getIndex(ni, nj, nk);
                        
                        // Add to queue if cavity and not visited
                        if (!grid_.occupied[nidx] && !visited[nidx]) {
                            queue.push(GridPoint(ni, nj, nk));
                            visited[nidx] = true;
                        }
                    }
                }
                
                // Calculate cluster properties
                if (count > 0) {
                    cluster.center = centerSum * (1.0 / count);
                    double spacingNm = gridSpacing_ * ANGSTROM_TO_NM;
                    cluster.volume = count * spacingNm * spacingNm * spacingNm;
                    
                    // Only keep clusters with significant volume
                    if (count >= 5) {  // At least 5 connected cells
                        clusters.push_back(cluster);
                    }
                }
            }
        }
    }
    
    // Sort clusters by volume (largest first)
    std::sort(clusters.begin(), clusters.end(), 
              [](const CavityCluster& a, const CavityCluster& b) {
                  return a.volume > b.volume;
              });
    
    // Log statistics
    stats_.clusterCount = static_cast<int>(clusters.size());
    if (!clusters.empty()) {
        stats_.largestClusterSize = static_cast<int>(clusters[0].positions.size());
        stats_.averageClusterSize = stats_.cavityPoints / std::max(1, stats_.clusterCount);
    }
    
    return clusters;
}

Vector3 CavityManager::selectFromCluster(const CavityCluster& cluster) {
    // Select random position from cluster
    if (cluster.positions.empty()) {
        return cluster.center;
    }
    
    int idx = utils::RandomUtils::uniformInt(0, static_cast<int>(cluster.positions.size()) - 1);
    return cluster.positions[idx];
}

std::vector<Vector3> CavityManager::getClusterCenters(const MCState& state, int maxClusters) {
    auto clusters = findCavityClustersFloodFill(state);
    
    if (clusters.empty()) {
        return std::vector<Vector3>();
    }
    
    std::vector<Vector3> centers;
    
    // Build cumulative distribution for size-weighted sampling
    std::vector<int> cumulative;
    cumulative.reserve(clusters.size() + 1);
    cumulative.push_back(0);
    
    for (const auto& cluster : clusters) {
        cumulative.push_back(cumulative.back() + static_cast<int>(cluster.positions.size()));
    }
    
    int totalPoints = cumulative.back();
    if (totalPoints == 0) {
        return centers;
    }
    
    // Select clusters with size-weighted probability
    std::unordered_set<int> selectedIndices;
    int numToSelect = std::min(maxClusters, static_cast<int>(clusters.size()));
    
    while (static_cast<int>(selectedIndices.size()) < numToSelect) {
        // Sample from weighted distribution
        int randomPoint = utils::RandomUtils::uniformInt(0, totalPoints - 1);
        
        // Find which cluster was selected
        int selectedCluster = 0;
        for (size_t i = 1; i < cumulative.size(); ++i) {
            if (randomPoint < cumulative[i]) {
                selectedCluster = i - 1;
                break;
            }
        }
        
        // Add to selected set (may already be selected, which is fine)
        selectedIndices.insert(selectedCluster);
    }
    
    // Extract centers from selected clusters
    for (int idx : selectedIndices) {
        centers.push_back(clusters[idx].center);
    }
    
    return centers;
}

// Color-based independent selection with size-weighted selection
std::vector<Vector3> CavityManager::selectIndependentCavities(
    const MCState& state, 
    double minSeparationNm,
    int maxPoints) {
    
    // Ensure cavities are up to date
    findCavities(state);
    
    // Calculate color grid size
    double spacingNm = gridSpacing_ * ANGSTROM_TO_NM;
    int colorSize = std::max(1, static_cast<int>(std::ceil(minSeparationNm / spacingNm)));
    
    // First, count cavities in each color class to enable weighted selection
    std::vector<int> colorClassSizes(colorSize * colorSize * colorSize, 0);
    std::vector<std::tuple<int, int, int>> colorIndices;
    
    // Count cavities in each color class
    for (int ci = 0; ci < colorSize; ++ci) {
        for (int cj = 0; cj < colorSize; ++cj) {
            for (int ck = 0; ck < colorSize; ++ck) {
                int count = 0;
                for (int i = ci; i < grid_.nx; i += colorSize) {
                    for (int j = cj; j < grid_.ny; j += colorSize) {
                        for (int k = ck; k < grid_.nz; k += colorSize) {
                            int idx = grid_.getIndex(i, j, k);
                            if (!grid_.occupied[idx]) {
                                count++;
                            }
                        }
                    }
                }
                int colorIdx = ci * colorSize * colorSize + cj * colorSize + ck;
                colorClassSizes[colorIdx] = count;
                if (count > 0) {
                    colorIndices.push_back(std::make_tuple(ci, cj, ck));
                }
            }
        }
    }
    
    // Select color class weighted by size (preserves uniform sampling over all cavities)
    if (colorIndices.empty()) {
        return std::vector<Vector3>();
    }
    
    // Build cumulative distribution for weighted sampling
    std::vector<int> cumulative;
    cumulative.reserve(colorIndices.size() + 1);
    cumulative.push_back(0);
    
    for (const auto& colorIdx : colorIndices) {
        int ci = std::get<0>(colorIdx);
        int cj = std::get<1>(colorIdx);
        int ck = std::get<2>(colorIdx);
        int idx = ci * colorSize * colorSize + cj * colorSize + ck;
        cumulative.push_back(cumulative.back() + colorClassSizes[idx]);
    }
    
    // Sample from weighted distribution
    int totalCavities = cumulative.back();
    if (totalCavities == 0) {
        return std::vector<Vector3>();
    }
    
    int randomCavity = utils::RandomUtils::uniformInt(0, totalCavities - 1);
    
    // Find which color class was selected
    int selectedClass = 0;
    for (size_t i = 1; i < cumulative.size(); ++i) {
        if (randomCavity < cumulative[i]) {
            selectedClass = i - 1;
            break;
        }
    }
    
    // Get the selected color indices
    int colorI = std::get<0>(colorIndices[selectedClass]);
    int colorJ = std::get<1>(colorIndices[selectedClass]);
    int colorK = std::get<2>(colorIndices[selectedClass]);
    
    std::vector<Vector3> selectedCavities;
    
    // Select all cavities in this color class
    for (int i = colorI; i < grid_.nx; i += colorSize) {
        for (int j = colorJ; j < grid_.ny; j += colorSize) {
            for (int k = colorK; k < grid_.nz; k += colorSize) {
                int idx = grid_.getIndex(i, j, k);
                
                if (!grid_.occupied[idx]) {
                    Vector3 posNm = gridToPosition(i, j, k);  // Already in nm
                    selectedCavities.push_back(posNm);
                    
                    if (static_cast<int>(selectedCavities.size()) >= maxPoints) {
                        return selectedCavities;
                    }
                }
            }
        }
    }
    
    return selectedCavities;
}

// Automatic mode selection based on system density
bool CavityManager::shouldUseClustering(const MCState& state) {
    // Use clustering if:
    // 1. Too many cavity points (>1000)
    // 2. Occupancy is moderate (30-70%)
    
    if (!cacheValid_) {
        findCavities(state);
    }
    
    bool tooManyCavities = stats_.cavityPoints > 1000;
    bool moderateOccupancy = stats_.occupancyRatio > 0.3 && stats_.occupancyRatio < 0.7;
    
    return tooManyCavities || moderateOccupancy;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc