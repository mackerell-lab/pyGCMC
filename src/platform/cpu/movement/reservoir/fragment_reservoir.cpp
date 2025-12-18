// Fragment Reservoir Stub Implementation
// Provides minimal implementation for Python binding tests

#include "../reservoir/fragment_reservoir.hpp"
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Statistics implementation
void FragmentReservoir::Statistics::print() const {
    std::cout << "=== FragmentReservoir Statistics ===" << std::endl;
    std::cout << "Global Statistics:" << std::endl;
    std::cout << "  Total Insertions: " << totalInsertions << std::endl;
    std::cout << "  Total Deletions: " << totalDeletions << std::endl;
    std::cout << "  Ghost Recycles: " << ghostRecycles << std::endl;
    std::cout << "  Memory Compactions: " << memoryCompactions << std::endl;
    std::cout << "  Peak Active Count: " << peakActiveCount << std::endl;
    std::cout << "  Average Active Count: " << averageActiveCount << std::endl;
    
    if (!activeCountByType.empty()) {
        std::cout << "\nPer-Type Statistics:" << std::endl;
        for (const auto& [typeId, count] : activeCountByType) {
            std::cout << "  Type " << typeId << ":" << std::endl;
            std::cout << "    Active: " << count;
            
            auto ghostIt = ghostCountByType.find(typeId);
            if (ghostIt != ghostCountByType.end()) {
                std::cout << ", Ghost: " << ghostIt->second;
            }
            
            auto lifetimeIt = averageLifetimeByType.find(typeId);
            if (lifetimeIt != averageLifetimeByType.end() && lifetimeIt->second > 0) {
                std::cout << ", Avg Lifetime: " << lifetimeIt->second;
            }
            
            auto acceptIt = acceptanceRateByType.find(typeId);
            if (acceptIt != acceptanceRateByType.end()) {
                std::cout << ", Accept Rate: " << (acceptIt->second * 100) << "%";
            }
            std::cout << std::endl;
        }
    }
    
    if (insertionTimeMs > 0 || deletionTimeMs > 0) {
        std::cout << "\nPerformance Metrics:" << std::endl;
        if (totalInsertions > 0) {
            std::cout << "  Avg Insertion Time: " << (insertionTimeMs / totalInsertions) << " ms" << std::endl;
        }
        if (totalDeletions > 0) {
            std::cout << "  Avg Deletion Time: " << (deletionTimeMs / totalDeletions) << " ms" << std::endl;
        }
        if (cacheHits + cacheMisses > 0) {
            double hitRate = static_cast<double>(cacheHits) / (cacheHits + cacheMisses);
            std::cout << "  Cache Hit Rate: " << (hitRate * 100) << "%" << std::endl;
        }
    }
}

void FragmentReservoir::Statistics::reset() {
    // Reset per-type statistics
    activeCountByType.clear();
    ghostCountByType.clear();
    averageLifetimeByType.clear();
    acceptanceRateByType.clear();
    
    // Reset global statistics
    totalInsertions = 0;
    totalDeletions = 0;
    ghostRecycles = 0;
    memoryCompactions = 0;
    averageGhostLifetime = 0.0;
    peakActiveCount = 0;
    averageActiveCount = 0.0;
    
    // Reset performance metrics
    insertionTimeMs = 0.0;
    deletionTimeMs = 0.0;
    queryTimeMs = 0.0;
    cacheHits = 0;
    cacheMisses = 0;
}

// FragmentReservoir stub implementation
FragmentReservoir::FragmentReservoir(std::shared_ptr<ActivePool> pool)
    : pool_(pool), config_(Config()), nextInstanceId_(0), nextTemplateId_(0), currentStep_(0) {
}

FragmentReservoir::FragmentReservoir(const Config& config, std::shared_ptr<ActivePool> pool)
    : pool_(pool), config_(config), nextInstanceId_(0), nextTemplateId_(0), currentStep_(0) {
}

FragmentReservoir::~FragmentReservoir() = default;

void FragmentReservoir::reserveInstanceIds(int startIndex) {
    nextInstanceId_ = std::max(nextInstanceId_, std::max(0, startIndex));
}

int FragmentReservoir::createInstanceWithId(int templateId,
                                            int instanceId,
                                            const Vector3& position,
                                            const Quaternion& orientation) {
    auto* tmpl = getTemplate(templateId);
    if (!tmpl) return -1;
    if (instanceId < 0) return -1;
    if (instances_.find(instanceId) != instances_.end()) return -1;

    // Ensure future auto-assigned IDs never collide with this fixed ID.
    nextInstanceId_ = std::max(nextInstanceId_, instanceId + 1);

    FragmentInstance instance;
    instance.templateId = templateId;
    instance.instanceId = instanceId;
    instance.position = position;
    instance.centerOfMass = position;
    instance.orientation = orientation;
    instance.isActive = true;
    instance.isGhost = false;
    instance.insertionTime = currentStep_;
    instance.residueIndex = instanceId;  // Stable mapping: instanceId == residueIndex

    instances_[instanceId] = instance;
    activeInstances_.insert(instanceId);
    templateInstances_[templateId].insert(instanceId);
    perTypeActiveCount_[templateId]++;

    return instanceId;
}

// Template management
int FragmentReservoir::addTemplate(const FragmentTemplate& tmpl) {
    int id = nextTemplateId_++;
    templates_[id] = tmpl;
    templateNameMap_[tmpl.name] = id;
    return id;
}

int FragmentReservoir::loadTemplate(const std::string& /*filename*/, const std::string& name, double chemicalPotential) {
    FragmentTemplate tmpl;
    tmpl.name = name;
    tmpl.chemicalPotential = chemicalPotential;
    return addTemplate(tmpl);
}

FragmentTemplate* FragmentReservoir::getTemplate(int templateId) {
    auto it = templates_.find(templateId);
    return (it != templates_.end()) ? &it->second : nullptr;
}

const FragmentTemplate* FragmentReservoir::getTemplate(int templateId) const {
    auto it = templates_.find(templateId);
    return (it != templates_.end()) ? &it->second : nullptr;
}

FragmentTemplate* FragmentReservoir::getTemplate(const std::string& name) {
    auto it = templateNameMap_.find(name);
    if (it != templateNameMap_.end()) {
        return getTemplate(it->second);
    }
    return nullptr;
}

const FragmentTemplate* FragmentReservoir::getTemplate(const std::string& name) const {
    auto it = templateNameMap_.find(name);
    if (it != templateNameMap_.end()) {
        return getTemplate(it->second);
    }
    return nullptr;
}

void FragmentReservoir::updateChemicalPotential(int templateId, double mu) {
    if (auto* tmpl = getTemplate(templateId)) {
        tmpl->chemicalPotential = mu;
    }
}

void FragmentReservoir::updateActivity(int templateId, double beta) {
    if (auto* tmpl = getTemplate(templateId)) {
        tmpl->updateActivity(beta);
    }
}

void FragmentReservoir::updateAllActivities(double beta) {
    for (auto& [id, tmpl] : templates_) {
        tmpl.updateActivity(beta);
    }
}

// Instance management
int FragmentReservoir::createInstance(int templateId, const Vector3& position, const Quaternion& orientation) {
    auto* tmpl = getTemplate(templateId);
    if (!tmpl) return -1;
    
    int instanceId = -1;
    
    // Priority 1: Try to recycle a ghost instance
    instanceId = recycleGhost(templateId);
    
    if (instanceId >= 0) {
        // Reuse existing ghost instance
        auto& instance = instances_[instanceId];
        instance.position = position;
        instance.centerOfMass = position;
        instance.orientation = orientation;
        instance.isActive = true;
        instance.isGhost = false;
        instance.insertionTime = currentStep_;
        // residueIndex will be set by pool if needed
        
        activeInstances_.insert(instanceId);
        ghostInstances_.erase(instanceId);
        stats_.ghostRecycles++;
    } else {
        // Allocate new instance
        instanceId = nextInstanceId_++;
        FragmentInstance instance;
        instance.templateId = templateId;
        instance.instanceId = instanceId;
        instance.position = position;
        instance.centerOfMass = position;
        instance.orientation = orientation;
        instance.isActive = true;
        instance.isGhost = false;
        instance.insertionTime = currentStep_;
        instance.residueIndex = -1;  // Initialize residueIndex
        
        instances_[instanceId] = instance;
        activeInstances_.insert(instanceId);
    }
    
    // Update template instances mapping
    templateInstances_[templateId].insert(instanceId);
    
    // Update per-type active count
    perTypeActiveCount_[templateId]++;
    
    // Update statistics
    stats_.totalInsertions++;
    
    // Update peak active count
    int currentActive = activeInstances_.size();
    if (currentActive > stats_.peakActiveCount) {
        stats_.peakActiveCount = currentActive;
    }
    
    // Update per-type statistics
    stats_.activeCountByType[templateId] = perTypeActiveCount_[templateId];
    
    // TODO: If pool_ is not null, call pool_->insertMolecule() to get residueIndex
    
    return instanceId;
}

int FragmentReservoir::createInstanceCBMC(int templateId, const std::vector<Vector3>& trialPositions,
                                         const std::vector<double>& /*trialEnergies*/) {
    if (trialPositions.empty()) return -1;
    return createInstance(templateId, trialPositions[0]);
}

bool FragmentReservoir::deleteInstance(int instanceId) {
    auto it = instances_.find(instanceId);
    if (it == instances_.end() || !it->second.isActive) return false;
    
    int templateId = it->second.templateId;
    
    // Mark as ghost
    it->second.isActive = false;
    it->second.isGhost = true;
    
    // Move from active to ghost sets
    activeInstances_.erase(instanceId);
    ghostInstances_.insert(instanceId);
    
    // Add to ghost pool for this template (FIFO)
    GhostRecord ghost;
    ghost.instanceId = instanceId;
    ghost.templateId = templateId;
    ghost.deletionStep = currentStep_;
    ghostPools_[templateId].push_back(ghost);
    
    // Update per-type active count
    if (perTypeActiveCount_[templateId] > 0) {
        perTypeActiveCount_[templateId]--;
    }
    
    // TODO: If pool_ is not null and residueIndex >= 0, call pool_->deleteResidue(residueIndex)
    int residueIndex = it->second.residueIndex;
    if (pool_ && residueIndex >= 0) {
        // pool_->deleteResidue(residueIndex); // Uncomment when ActivePool is ready
        it->second.residueIndex = -1;
    }
    
    // Update statistics
    stats_.totalDeletions++;
    
    // Update per-type statistics
    stats_.activeCountByType[templateId] = perTypeActiveCount_[templateId];
    stats_.ghostCountByType[templateId] = ghostPools_[templateId].size();
    
    return true;
}

bool FragmentReservoir::restoreInstance(int instanceId, const Vector3& position, const Quaternion& orientation) {
    auto it = instances_.find(instanceId);
    if (it == instances_.end() || it->second.isActive) return false;
    
    int templateId = it->second.templateId;
    
    // Reset instance state
    it->second.isActive = true;
    it->second.isGhost = false;
    it->second.position = position;
    it->second.centerOfMass = position;
    it->second.orientation = orientation;
    it->second.insertionTime = currentStep_;
    
    // Move from ghost to active sets
    ghostInstances_.erase(instanceId);
    activeInstances_.insert(instanceId);
    
    // Remove from ghost pool if present
    auto poolIt = ghostPools_.find(templateId);
    if (poolIt != ghostPools_.end()) {
        auto& ghosts = poolIt->second;
        ghosts.erase(
            std::remove_if(ghosts.begin(), ghosts.end(),
                          [instanceId](const GhostRecord& g) { return g.instanceId == instanceId; }),
            ghosts.end()
        );
    }
    
    // Update per-type active count
    perTypeActiveCount_[templateId]++;
    
    // TODO: If pool_ is not null, call pool_->insertMolecule() to reallocate residueIndex
    
    return true;
}

bool FragmentReservoir::purgeInstance(int instanceId) {
    auto it = instances_.find(instanceId);
    if (it == instances_.end()) return false;
    
    activeInstances_.erase(instanceId);
    ghostInstances_.erase(instanceId);
    instances_.erase(it);
    
    return true;
}

std::vector<int> FragmentReservoir::createMultipleInstances(int templateId, const std::vector<Vector3>& positions) {
    std::vector<int> ids;
    for (const auto& pos : positions) {
        int id = createInstance(templateId, pos);
        if (id >= 0) ids.push_back(id);
    }
    return ids;
}

int FragmentReservoir::deleteMultipleInstances(const std::vector<int>& instanceIds) {
    int count = 0;
    for (int id : instanceIds) {
        if (deleteInstance(id)) count++;
    }
    return count;
}

// Ghost management
int FragmentReservoir::recycleGhost(int templateId) {
    // Check if there are ghosts for this template
    auto it = ghostPools_.find(templateId);
    if (it == ghostPools_.end() || it->second.empty()) {
        return -1;  // No ghosts available
    }
    
    // Get the oldest ghost (FIFO)
    GhostRecord ghost = it->second.front();
    it->second.pop_front();
    
    // Remove from ghostInstances_ set
    ghostInstances_.erase(ghost.instanceId);
    
    // Return the instance ID for reuse
    return ghost.instanceId;
}

int FragmentReservoir::purgeGhosts(int maxToKeep) {
    // Use config maxGhosts if maxToKeep is -1
    if (maxToKeep < 0) {
        maxToKeep = config_.maxGhosts;
    }
    
    int totalPurged = 0;
    
    // Process each template's ghost pool
    for (auto& [templateId, ghosts] : ghostPools_) {
        // Keep only the newest maxToKeep ghosts
        while (static_cast<int>(ghosts.size()) > maxToKeep) {
            // Remove the oldest ghost (from front)
            GhostRecord ghost = ghosts.front();
            ghosts.pop_front();
            
            // Completely purge this instance
            purgeInstance(ghost.instanceId);
            totalPurged++;
        }
    }
    
    return totalPurged;
}

int FragmentReservoir::getGhostCount(int templateId) const {
    if (templateId < 0) {
        // Return total ghost count across all types
        int totalCount = 0;
        for (const auto& [tid, ghosts] : ghostPools_) {
            totalCount += ghosts.size();
        }
        return totalCount;
    }
    
    // Return ghost count for specific template
    auto it = ghostPools_.find(templateId);
    if (it != ghostPools_.end()) {
        return it->second.size();
    }
    return 0;
}

// Query methods
FragmentInstance* FragmentReservoir::getInstance(int instanceId) {
    auto it = instances_.find(instanceId);
    return (it != instances_.end()) ? &it->second : nullptr;
}

const FragmentInstance* FragmentReservoir::getInstance(int instanceId) const {
    auto it = instances_.find(instanceId);
    return (it != instances_.end()) ? &it->second : nullptr;
}

FragmentInstance* FragmentReservoir::getInstanceByResidueIndex(int residueIdx) {
    for (auto& [id, inst] : instances_) {
        if (inst.residueIndex == residueIdx && inst.isActive) {
            return &inst;
        }
    }
    return nullptr;
}

std::vector<int> FragmentReservoir::getActiveInstances(int templateId) const {
    std::vector<int> result;
    if (templateId < 0) {
        result.assign(activeInstances_.begin(), activeInstances_.end());
    } else {
        for (int id : activeInstances_) {
            auto it = instances_.find(id);
            if (it != instances_.end() && it->second.templateId == templateId) {
                result.push_back(id);
            }
        }
    }
    return result;
}

int FragmentReservoir::getActiveCount(int templateId) const {
    if (templateId < 0) return activeInstances_.size();
    
    int count = 0;
    for (int id : activeInstances_) {
        auto it = instances_.find(id);
        if (it != instances_.end() && it->second.templateId == templateId) {
            count++;
        }
    }
    return count;
}

std::vector<int> FragmentReservoir::findInstancesInSphere(const Vector3& center, double radius) const {
    std::vector<int> result;
    double radiusSq = radius * radius;
    
    for (int id : activeInstances_) {
        auto it = instances_.find(id);
        if (it != instances_.end()) {
            const auto& pos = it->second.position;
            double distSq = (pos.x - center.x) * (pos.x - center.x) +
                          (pos.y - center.y) * (pos.y - center.y) +
                          (pos.z - center.z) * (pos.z - center.z);
            if (distSq <= radiusSq) {
                result.push_back(id);
            }
        }
    }
    return result;
}

std::vector<int> FragmentReservoir::findInstancesInBox(const Vector3& min, const Vector3& max) const {
    std::vector<int> result;
    
    for (int id : activeInstances_) {
        auto it = instances_.find(id);
        if (it != instances_.end()) {
            const auto& pos = it->second.position;
            if (pos.x >= min.x && pos.x <= max.x &&
                pos.y >= min.y && pos.y <= max.y &&
                pos.z >= min.z && pos.z <= max.z) {
                result.push_back(id);
            }
        }
    }
    return result;
}

std::vector<int> FragmentReservoir::getInstancesByEnergy(double minE, double maxE) const {
    std::vector<int> result;
    
    for (int id : activeInstances_) {
        auto it = instances_.find(id);
        if (it != instances_.end()) {
            double energy = it->second.energy_total;
            if (energy >= minE && energy <= maxE) {
                result.push_back(id);
            }
        }
    }
    return result;
}

std::vector<int> FragmentReservoir::getInstancesByLifetime(double minTime, double maxTime, double currentStep) const {
    std::vector<int> result;
    
    for (int id : activeInstances_) {
        auto it = instances_.find(id);
        if (it != instances_.end()) {
            double lifetime = it->second.getLifetime(currentStep);
            if (lifetime >= minTime && lifetime <= maxTime) {
                result.push_back(id);
            }
        }
    }
    return result;
}

std::vector<int> FragmentReservoir::getGhostIndices() const {
    return std::vector<int>(ghostInstances_.begin(), ghostInstances_.end());
}

// Note: getStatistics returns by reference, already defined inline in header

void FragmentReservoir::resetStatistics() {
    stats_.reset();
}

void FragmentReservoir::printStatistics() const {
    stats_.print();
}

// Memory management
void FragmentReservoir::compact() {
    // Check if compaction is needed based on config
    if (!config_.autoCompact) return;
    
    double fragmentation = getFragmentation();
    if (fragmentation < config_.compactThreshold) return;
    
    // If pool is available, delegate compaction to it
    if (pool_) {
        // pool_->compact(); // Uncomment when ActivePool::compact is ready
        stats_.memoryCompactions++;
    }
    
    // Clean up excess ghosts if needed
    if (config_.maxGhosts > 0) {
        purgeGhosts(config_.maxGhosts);
    }
}

double FragmentReservoir::getFragmentation() const {
    if (!pool_) return 0.0;
    
    // Simple fragmentation metric: ratio of ghosts to total instances
    int totalInstances = instances_.size();
    if (totalInstances == 0) return 0.0;
    
    int ghostCount = ghostInstances_.size();
    return static_cast<double>(ghostCount) / totalInstances;
    
    // TODO: When ActivePool provides getFragmentation(), use:
    // return pool_->getFragmentation();
}

// Additional stub methods required by GCMCEngine
void FragmentReservoir::updatePosition(int instanceId, const Vector3& newPos) {
    auto it = instances_.find(instanceId);
    if (it != instances_.end()) {
        it->second.position = newPos;
        it->second.centerOfMass = newPos;
    }
}

void FragmentReservoir::updateOrientation(int instanceId, const Quaternion& newOrient) {
    auto it = instances_.find(instanceId);
    if (it != instances_.end()) {
        it->second.orientation = newOrient;
    }
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
