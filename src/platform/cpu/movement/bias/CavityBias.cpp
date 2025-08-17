#include "CavityBias.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../common/MovementUtils.hpp"
#include <cmath>
#include <algorithm>
#include <queue>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

// Constants for unit conversion
static constexpr double ANGSTROM_TO_NM = 0.1;
static constexpr double NM_TO_ANGSTROM = 10.0;

CavityManager::CavityManager(double gridSpacing, double probeRadius)
    : gridSpacing_(gridSpacing),      // in nm
      probeRadius_(probeRadius),      // in nm
      cacheValid_(false) {
    resetStatistics();
}

CavityManager::~CavityManager() = default;

std::vector<Vector3> CavityManager::findCavities(const MCState& state) {
    // Use cache if valid
    if (cacheValid_ && !cavityCache_.empty()) {
        stats_.cacheHits++;
        return cavityCache_;
    }
    
    stats_.cacheMisses++;
    
    // Get box dimensions (in nm from MCState)
    Vector3 boxSize(state.info.box[0], state.info.box[1], state.info.box[2]);
    
    // Initialize grid
    initializeGrid(boxSize);
    
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
                    Vector3 pos = gridToPosition(i, j, k);
                    // Check if it's truly a cavity (not just unoccupied)
                    if (checkCavity(pos, state)) {
                        cavityCache_.push_back(pos);
                    }
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
    // Convert position from nm to grid coordinates
    Vector3 posAngstrom = position * NM_TO_ANGSTROM;
    
    int i = static_cast<int>((posAngstrom.x - grid_.origin.x) / grid_.spacing.x);
    int j = static_cast<int>((posAngstrom.y - grid_.origin.y) / grid_.spacing.y);
    int k = static_cast<int>((posAngstrom.z - grid_.origin.z) / grid_.spacing.z);
    
    if (!grid_.isValid(i, j, k)) {
        return false;
    }
    
    return !grid_.occupied[grid_.getIndex(i, j, k)];
}

void CavityManager::invalidateCache() {
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
    double clusterThreshold = gridSpacing_ * 1.5 * ANGSTROM_TO_NM;  // Convert to nm
    
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
            newCluster.volume = std::pow(gridSpacing_ * ANGSTROM_TO_NM, 3);  // Initial volume
            clusters.push_back(newCluster);
        }
    }
    
    // Calculate cluster volumes
    for (auto& cluster : clusters) {
        cluster.volume = cluster.positions.size() * std::pow(gridSpacing_ * ANGSTROM_TO_NM, 3);
    }
    
    return clusters;
}

// Private helper functions

void CavityManager::initializeGrid(const Vector3& boxSize) {
    // Box size is in nm, convert grid spacing to nm for calculation
    double spacingNm = gridSpacing_ * ANGSTROM_TO_NM;
    
    // Calculate grid dimensions
    grid_.nx = std::max(1, static_cast<int>(std::ceil(boxSize.x / spacingNm)));
    grid_.ny = std::max(1, static_cast<int>(std::ceil(boxSize.y / spacingNm)));
    grid_.nz = std::max(1, static_cast<int>(std::ceil(boxSize.z / spacingNm)));
    
    // Set grid properties (store in Angstroms for consistency)
    grid_.origin = Vector3(0.0, 0.0, 0.0);
    grid_.spacing = Vector3(gridSpacing_, gridSpacing_, gridSpacing_);
    grid_.boxSize = boxSize * NM_TO_ANGSTROM;  // Convert to Angstroms
    
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
                // Position is in nm, convert to Angstroms
                Vector3 posAngstrom(atom.x * NM_TO_ANGSTROM, 
                                   atom.y * NM_TO_ANGSTROM, 
                                   atom.z * NM_TO_ANGSTROM);
                
                // Use default VDW radius if not specified (1.5 Å for simplicity)
                double radius = 1.5 + probeRadius_;  // in Angstroms
                markOccupiedRegion(posAngstrom, radius);
            }
        }
    }
}

void CavityManager::markOccupiedRegion(const Vector3& center, double radius) {
    // Center and radius are in Angstroms
    int iMin = std::max(0, static_cast<int>((center.x - radius - grid_.origin.x) / grid_.spacing.x));
    int iMax = std::min(grid_.nx - 1, static_cast<int>((center.x + radius - grid_.origin.x) / grid_.spacing.x));
    int jMin = std::max(0, static_cast<int>((center.y - radius - grid_.origin.y) / grid_.spacing.y));
    int jMax = std::min(grid_.ny - 1, static_cast<int>((center.y + radius - grid_.origin.y) / grid_.spacing.y));
    int kMin = std::max(0, static_cast<int>((center.z - radius - grid_.origin.z) / grid_.spacing.z));
    int kMax = std::min(grid_.nz - 1, static_cast<int>((center.z + radius - grid_.origin.z) / grid_.spacing.z));
    
    double radiusSq = radius * radius;
    
    for (int i = iMin; i <= iMax; ++i) {
        for (int j = jMin; j <= jMax; ++j) {
            for (int k = kMin; k <= kMax; ++k) {
                Vector3 gridPos = gridToPosition(i, j, k);  // Already in nm
                Vector3 diff = gridPos - center;
                double distSq = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
                
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
    double minDist = probeRadius_ * ANGSTROM_TO_NM;  // Convert to nm
    double minDistSq = minDist * minDist;
    
    for (int resIdx = 0; resIdx < state.activeResidueCount; ++resIdx) {
        const MCResidue& residue = state.residues[resIdx];
        if (!residue.active) continue;
        
        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                const MCAtom& atom = state.atoms[atomIdx];
                // Calculate distance squared with PBC
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
    
    // Fall back to random position
    stats_.randomInsertions++;
    return selectRandomPosition(Vector3(state.info.box[0], state.info.box[1], state.info.box[2]));
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

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc