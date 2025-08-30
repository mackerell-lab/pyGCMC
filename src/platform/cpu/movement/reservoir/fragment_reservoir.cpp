// Fragment Reservoir System Implementation
#include "fragment_reservoir.hpp"
#include "../../../../model/montecarlo/MCStructures.hpp"
#include "../pool/ActivePool.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <numeric>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCAtom;
using model::montecarlo::MCResidue;

// ============================================================================
// Constructor/Destructor
// ============================================================================

FragmentReservoir::FragmentReservoir(std::shared_ptr<ActivePool> pool)
    : activePool_(pool), ownsPool_(false) {
    if (!activePool_) {
        // Create our own pool if none provided
        activePool_ = std::make_shared<ActivePool>(100000, 30000);
        ownsPool_ = true;
    }
    instances_.reserve(config_.maxInstances);
}

FragmentReservoir::FragmentReservoir(const Config& config, 
                                   std::shared_ptr<ActivePool> pool)
    : activePool_(pool), ownsPool_(false), config_(config) {
    if (!activePool_) {
        activePool_ = std::make_shared<ActivePool>(100000, 30000);
        ownsPool_ = true;
    }
    instances_.reserve(config_.maxInstances);
}

FragmentReservoir::~FragmentReservoir() {
    // Clean up any remaining ghosts
    if (ownsPool_) {
        purgeGhosts(0);
    }
}

// ============================================================================
// Template Management
// ============================================================================

int FragmentReservoir::addTemplate(const FragmentTemplate& tmpl) {
    int id = static_cast<int>(templates_.size());
    templates_.push_back(tmpl);
    templateNameMap_[tmpl.name] = id;
    typeToTemplateMap_[tmpl.typeId] = id;
    
    // Initialize statistics for this type
    stats_.activeCountByType[id] = 0;
    stats_.ghostCountByType[id] = 0;
    stats_.averageLifetimeByType[id] = 0.0;
    stats_.acceptanceRateByType[id] = 0.0;
    
    return id;
}

int FragmentReservoir::loadTemplate(const std::string& filename, 
                                   const std::string& name,
                                   double chemicalPotential) {
    FragmentTemplate tmpl;
    tmpl.name = name;
    tmpl.chemicalPotential = chemicalPotential;
    
    // TODO: Implement PDB/MOL2 parsing from filename
    (void)filename;  // Suppress unused parameter warning until implemented
    // For now, create a simple water molecule as example
    if (name == "WAT" || name == "TIP3") {
        // Water molecule
        tmpl.molecularWeight = 18.015;
        tmpl.radius = 1.4;  // Å
        
        // Add atoms (O, H, H)
        MCAtom oxygen;
        oxygen.x = 0.0; oxygen.y = 0.0; oxygen.z = 0.0;
        oxygen.charge = -0.834;
        oxygen.type = 0;  // Oxygen type
        tmpl.atoms.push_back(oxygen);
        
        MCAtom hydrogen1;
        hydrogen1.x = 0.9572; hydrogen1.y = 0.0; hydrogen1.z = 0.0;
        hydrogen1.charge = 0.417;
        hydrogen1.type = 1;  // Hydrogen type
        tmpl.atoms.push_back(hydrogen1);
        
        MCAtom hydrogen2;
        hydrogen2.x = -0.2399; hydrogen2.y = 0.9266; hydrogen2.z = 0.0;
        hydrogen2.charge = 0.417;
        hydrogen2.type = 1;
        tmpl.atoms.push_back(hydrogen2);
        
        // Add bonds
        tmpl.bonds.push_back({0, 1, 0.9572});
        tmpl.bonds.push_back({0, 2, 0.9572});
        
        // Add angle
        tmpl.angles.push_back({1, 0, 2, 104.52 * M_PI / 180.0});
    }
    
    tmpl.typeId = static_cast<int>(templates_.size());
    return addTemplate(tmpl);
}

FragmentTemplate* FragmentReservoir::getTemplate(int templateId) {
    if (isValidTemplateId(templateId)) {
        return &templates_[templateId];
    }
    return nullptr;
}

const FragmentTemplate* FragmentReservoir::getTemplate(int templateId) const {
    if (isValidTemplateId(templateId)) {
        return &templates_[templateId];
    }
    return nullptr;
}

FragmentTemplate* FragmentReservoir::getTemplate(const std::string& name) {
    auto it = templateNameMap_.find(name);
    if (it != templateNameMap_.end()) {
        return &templates_[it->second];
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
    for (auto& tmpl : templates_) {
        tmpl.updateActivity(beta);
    }
}

// ============================================================================
// Instance Creation
// ============================================================================

int FragmentReservoir::createInstance(int templateId, 
                                     const Vector3& position,
                                     const Quaternion& orientation) {
    if (!isValidTemplateId(templateId)) {
        return -1;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Try to recycle a ghost first
    int instanceSlot = -1;
    if (config_.ghostRecycleRatio > 0 && !ghostIndices_.empty()) {
        double ghostRatio = static_cast<double>(ghostIndices_.size()) / 
                          (instances_.size() + 1);
        if (ghostRatio > config_.ghostRecycleRatio) {
            instanceSlot = recycleGhost(templateId);
        }
    }
    
    // Allocate new slot if no ghost was recycled
    if (instanceSlot < 0) {
        instanceSlot = allocateInstanceSlot();
        if (instanceSlot < 0) {
            return -1;  // No space available
        }
    }
    
    // Get template
    const FragmentTemplate& tmpl = templates_[templateId];
    FragmentInstance& instance = instances_[instanceSlot];
    
    // Initialize instance
    instance.templateId = templateId;
    instance.instanceId = nextInstanceId_++;
    instance.isActive = true;
    instance.isGhost = false;
    instance.centerOfMass = position;
    instance.orientation = orientation;
    instance.insertionTime = currentMCStep_;
    instance.lastMoveTime = currentMCStep_;
    instance.energy_vdw = 0.0;
    instance.energy_elec = 0.0;
    instance.energy_total = 0.0;
    
    // Transform atoms to world coordinates
    std::vector<MCAtom> transformedAtoms = transformAtoms(tmpl.atoms, position, orientation);
    
    // Insert into ActivePool
    if (activePool_) {
        instance.residueIndex = activePool_->insertMolecule(transformedAtoms, tmpl.typeId);
        if (instance.residueIndex < 0) {
            // Failed to insert into pool
            freeInstanceSlot(instanceSlot);
            return -1;
        }
    }
    
    // Update indices
    activeByType_[templateId].push_back(instanceSlot);
    
    // Update statistics
    stats_.totalInsertions++;
    stats_.activeCountByType[templateId]++;
    if (stats_.activeCountByType[templateId] > stats_.peakActiveCount) {
        stats_.peakActiveCount = stats_.activeCountByType[templateId];
    }
    
    const_cast<FragmentTemplate&>(tmpl).totalInsertions++;
    
    auto end = std::chrono::high_resolution_clock::now();
    stats_.insertionTimeMs += std::chrono::duration<double, std::milli>(end - start).count();
    
    return instanceSlot;
}

int FragmentReservoir::createInstanceCBMC(int templateId,
                                         const std::vector<Vector3>& trialPositions,
                                         const std::vector<double>& trialEnergies) {
    if (!isValidTemplateId(templateId) || trialPositions.empty()) {
        return -1;
    }
    
    // Calculate Rosenbluth weight
    double W_new = 0.0;
    for (double energy : trialEnergies) {
        W_new += std::exp(-energy);  // Assuming β=1 for now
    }
    
    // Select configuration based on Boltzmann weights
    std::vector<double> probabilities;
    for (double energy : trialEnergies) {
        probabilities.push_back(std::exp(-energy) / W_new);
    }
    
    // Random selection based on probabilities
    double r = static_cast<double>(rand()) / RAND_MAX;
    double cumsum = 0.0;
    int selected = 0;
    for (size_t i = 0; i < probabilities.size(); ++i) {
        cumsum += probabilities[i];
        if (r < cumsum) {
            selected = i;
            break;
        }
    }
    
    // Create instance at selected position
    return createInstance(templateId, trialPositions[selected]);
}

// ============================================================================
// Instance Deletion
// ============================================================================

bool FragmentReservoir::deleteInstance(int instanceId) {
    if (!isValidInstanceId(instanceId)) {
        return false;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    FragmentInstance& instance = instances_[instanceId];
    if (instance.isGhost) {
        return false;  // Already a ghost
    }
    
    int templateId = instance.templateId;
    
    // Convert to ghost
    convertToGhost(instanceId);
    
    // Update statistics
    stats_.totalDeletions++;
    stats_.activeCountByType[templateId]--;
    stats_.ghostCountByType[templateId]++;
    
    // Update lifetime
    double lifetime = currentMCStep_ - instance.insertionTime;
    updateAverageLifetime(templateId, lifetime);
    
    auto end = std::chrono::high_resolution_clock::now();
    stats_.deletionTimeMs += std::chrono::duration<double, std::milli>(end - start).count();
    
    // Auto-purge if too many ghosts
    if (static_cast<int>(ghostIndices_.size()) > config_.maxGhosts) {
        purgeGhosts(config_.maxGhosts);
    }
    
    return true;
}

bool FragmentReservoir::purgeInstance(int instanceId) {
    if (!isValidInstanceId(instanceId)) {
        return false;
    }
    
    FragmentInstance& instance = instances_[instanceId];
    
    // Remove from ActivePool
    if (activePool_ && instance.residueIndex >= 0) {
        activePool_->deleteResidue(instance.residueIndex);
    }
    
    // Free the slot
    freeInstanceSlot(instanceId);
    
    // Update ghost statistics if it was a ghost
    if (instance.isGhost) {
        ghostIndices_.erase(instanceId);
        stats_.ghostCountByType[instance.templateId]--;
    }
    
    return true;
}

// ============================================================================
// Ghost Management
// ============================================================================

int FragmentReservoir::recycleGhost(int templateId) {
    int selected = selectGhostForRecycling(templateId);
    if (selected >= 0) {
        recycleGhostSlot(selected);
        stats_.ghostRecycles++;
        return selected;
    }
    return -1;
}

int FragmentReservoir::purgeGhosts(int maxToKeep) {
    if (maxToKeep < 0) {
        maxToKeep = config_.maxGhosts;
    }
    
    int purged = 0;
    while (static_cast<int>(ghostIndices_.size()) > maxToKeep && !ghostIndices_.empty()) {
        int idx = *ghostIndices_.begin();
        if (purgeInstance(idx)) {
            purged++;
        }
    }
    
    return purged;
}

int FragmentReservoir::getGhostCount(int templateId) const {
    if (templateId < 0) {
        return static_cast<int>(ghostIndices_.size());
    }
    
    auto it = stats_.ghostCountByType.find(templateId);
    return (it != stats_.ghostCountByType.end()) ? it->second : 0;
}

// ============================================================================
// Query Operations
// ============================================================================

std::vector<int> FragmentReservoir::getActiveInstances(int templateId) const {
    std::vector<int> result;
    
    if (templateId >= 0) {
        auto it = activeByType_.find(templateId);
        if (it != activeByType_.end()) {
            result = it->second;
        }
    } else {
        for (const auto& [type, indices] : activeByType_) {
            result.insert(result.end(), indices.begin(), indices.end());
        }
    }
    
    return result;
}

int FragmentReservoir::getActiveCount(int templateId) const {
    if (templateId >= 0) {
        auto it = stats_.activeCountByType.find(templateId);
        return (it != stats_.activeCountByType.end()) ? it->second : 0;
    }
    
    int total = 0;
    for (const auto& [type, count] : stats_.activeCountByType) {
        total += count;
    }
    return total;
}

FragmentInstance* FragmentReservoir::getInstance(int instanceId) {
    if (isValidInstanceId(instanceId)) {
        return &instances_[instanceId];
    }
    return nullptr;
}

const FragmentInstance* FragmentReservoir::getInstance(int instanceId) const {
    if (isValidInstanceId(instanceId)) {
        return &instances_[instanceId];
    }
    return nullptr;
}

std::vector<int> FragmentReservoir::findInstancesInSphere(const Vector3& center, 
                                                         double radius) const {
    std::vector<int> result;
    double r2 = radius * radius;
    
    for (size_t i = 0; i < instances_.size(); ++i) {
        const auto& inst = instances_[i];
        if (!inst.isActive || inst.isGhost) continue;
        
        Vector3 diff = inst.centerOfMass - center;
        double dist2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
        if (dist2 <= r2) {
            result.push_back(i);
        }
    }
    
    return result;
}

// ============================================================================
// Energy Management
// ============================================================================

void FragmentReservoir::updateEnergy(int instanceId, double vdw, double elec) {
    if (auto* inst = getInstance(instanceId)) {
        inst->energy_vdw = vdw;
        inst->energy_elec = elec;
        inst->energy_total = vdw + elec;
        inst->lastEnergyUpdate = currentMCStep_;
    }
}

void FragmentReservoir::updateTotalEnergy(int instanceId, double total) {
    if (auto* inst = getInstance(instanceId)) {
        inst->energy_total = total;
        inst->lastEnergyUpdate = currentMCStep_;
    }
}

// ============================================================================
// Movement Tracking
// ============================================================================

void FragmentReservoir::recordMoveAttempt(int instanceId, bool accepted, 
                                         double currentStep) {
    if (auto* inst = getInstance(instanceId)) {
        inst->moveAttempts++;
        if (accepted) {
            inst->acceptedMoves++;
            inst->lastMoveTime = currentStep;
        }
        
        // Update template statistics
        int templateId = inst->templateId;
        double rate = inst->getAcceptanceRate();
        updateAcceptanceRate(templateId, rate);
    }
    currentMCStep_ = currentStep;
}

void FragmentReservoir::updatePosition(int instanceId, const Vector3& newPos) {
    if (auto* inst = getInstance(instanceId)) {
        inst->centerOfMass = newPos;
    }
}

void FragmentReservoir::updateOrientation(int instanceId, const Quaternion& newOrient) {
    if (auto* inst = getInstance(instanceId)) {
        inst->orientation = newOrient;
    }
}

// ============================================================================
// Statistics
// ============================================================================

void FragmentReservoir::Statistics::print() const {
    std::cout << "\n=== Fragment Reservoir Statistics ===" << std::endl;
    std::cout << "Total insertions: " << totalInsertions << std::endl;
    std::cout << "Total deletions: " << totalDeletions << std::endl;
    std::cout << "Ghost recycles: " << ghostRecycles << std::endl;
    std::cout << "Memory compactions: " << memoryCompactions << std::endl;
    std::cout << "Peak active count: " << peakActiveCount << std::endl;
    
    std::cout << "\nPer-type statistics:" << std::endl;
    for (const auto& [type, count] : activeCountByType) {
        std::cout << "  Type " << type << ": " << count << " active";
        auto it = ghostCountByType.find(type);
        if (it != ghostCountByType.end() && it->second > 0) {
            std::cout << ", " << it->second << " ghosts";
        }
        auto it2 = averageLifetimeByType.find(type);
        if (it2 != averageLifetimeByType.end()) {
            std::cout << ", avg lifetime: " << it2->second;
        }
        std::cout << std::endl;
    }
    
    std::cout << "\nPerformance:" << std::endl;
    std::cout << "  Avg insertion time: " << insertionTimeMs / totalInsertions << " ms" << std::endl;
    std::cout << "  Avg deletion time: " << deletionTimeMs / totalDeletions << " ms" << std::endl;
    if (cacheHits + cacheMisses > 0) {
        double hitRate = static_cast<double>(cacheHits) / (cacheHits + cacheMisses);
        std::cout << "  Cache hit rate: " << hitRate * 100 << "%" << std::endl;
    }
}

void FragmentReservoir::Statistics::reset() {
    activeCountByType.clear();
    ghostCountByType.clear();
    averageLifetimeByType.clear();
    acceptanceRateByType.clear();
    totalInsertions = 0;
    totalDeletions = 0;
    ghostRecycles = 0;
    memoryCompactions = 0;
    averageGhostLifetime = 0.0;
    peakActiveCount = 0;
    averageActiveCount = 0.0;
    insertionTimeMs = 0.0;
    deletionTimeMs = 0.0;
    queryTimeMs = 0.0;
    cacheHits = 0;
    cacheMisses = 0;
}

void FragmentReservoir::printStatistics() const {
    stats_.print();
}

// ============================================================================
// Private Helper Functions
// ============================================================================

int FragmentReservoir::allocateInstanceSlot() {
    // Try to reuse a free slot
    if (!freeInstanceSlots_.empty()) {
        int slot = freeInstanceSlots_.front();
        freeInstanceSlots_.pop();
        return slot;
    }
    
    // Allocate new slot
    if (instances_.size() < static_cast<size_t>(config_.maxInstances)) {
        int slot = instances_.size();
        instances_.emplace_back();
        return slot;
    }
    
    return -1;  // No space available
}

void FragmentReservoir::freeInstanceSlot(int slot) {
    if (slot >= 0 && slot < static_cast<int>(instances_.size())) {
        instances_[slot] = FragmentInstance();  // Reset
        instances_[slot].instanceId = -1;  // Mark as free
        freeInstanceSlots_.push(slot);
    }
}

std::vector<MCAtom> FragmentReservoir::transformAtoms(const std::vector<MCAtom>& atoms,
                                                      const Vector3& position,
                                                      const Quaternion& orientation) const {
    std::vector<MCAtom> transformed = atoms;
    
    // Apply rotation and translation
    for (auto& atom : transformed) {
        Vector3 pos(atom.x, atom.y, atom.z);
        
        // Quaternion rotation: v' = q * v * q^-1
        // Using the formula for rotating a vector by a quaternion
        double qw = orientation.w, qx = orientation.x, qy = orientation.y, qz = orientation.z;
        
        // Cross product terms
        double tx = 2.0 * (qy * pos.z - qz * pos.y);
        double ty = 2.0 * (qz * pos.x - qx * pos.z);
        double tz = 2.0 * (qx * pos.y - qy * pos.x);
        
        // Apply rotation
        pos.x += qw * tx + qy * tz - qz * ty;
        pos.y += qw * ty + qz * tx - qx * tz;
        pos.z += qw * tz + qx * ty - qy * tx;
        
        // Apply translation
        pos = pos + position;
        
        atom.x = pos.x;
        atom.y = pos.y;
        atom.z = pos.z;
    }
    
    return transformed;
}

void FragmentReservoir::convertToGhost(int instanceId) {
    FragmentInstance& instance = instances_[instanceId];
    instance.isGhost = true;
    instance.isActive = false;
    
    // Remove from active list
    auto& activeList = activeByType_[instance.templateId];
    activeList.erase(
        std::remove(activeList.begin(), activeList.end(), instanceId),
        activeList.end()
    );
    
    // Add to ghost set
    ghostIndices_.insert(instanceId);
    
    // Mark as inactive in ActivePool
    if (activePool_ && instance.residueIndex >= 0) {
        if (auto* meta = activePool_->getResidueMetadata(instance.residueIndex)) {
            meta->active = false;
        }
    }
}

void FragmentReservoir::recycleGhostSlot(int instanceId) {
    FragmentInstance& instance = instances_[instanceId];
    
    // Note: The instance will be fully reinitialized by createInstance
    // Here we just update the ghost status
    
    // Remove from ghost set
    ghostIndices_.erase(instanceId);
    stats_.ghostCountByType[instance.templateId]--;
    
    // Note: Don't set isActive/isGhost here - let createInstance handle it
    // Don't add to activeByType_ here - createInstance will do it
}

int FragmentReservoir::selectGhostForRecycling(int templateId) {
    // Prefer ghosts of the same type
    for (int idx : ghostIndices_) {
        if (instances_[idx].templateId == templateId) {
            return idx;
        }
    }
    
    // Use any ghost if no type match
    if (!ghostIndices_.empty()) {
        return *ghostIndices_.begin();
    }
    
    return -1;
}

void FragmentReservoir::updateAverageLifetime(int templateId, double lifetime) {
    auto& avg = stats_.averageLifetimeByType[templateId];
    const auto& tmpl = templates_[templateId];
    int n = tmpl.totalInsertions;
    
    if (n > 0) {
        avg = (avg * (n - 1) + lifetime) / n;
    } else {
        avg = lifetime;
    }
}

void FragmentReservoir::updateAcceptanceRate(int templateId, double rate) {
    auto& avg = stats_.acceptanceRateByType[templateId];
    static std::map<int, int> updateCounts;
    updateCounts[templateId]++;
    int n = updateCounts[templateId];
    
    if (n > 0) {
        avg = (avg * (n - 1) + rate) / n;
    } else {
        avg = rate;
    }
}

// ============================================================================
// Memory Management
// ============================================================================

void FragmentReservoir::compact() {
    if (activePool_) {
        int compacted = activePool_->compact(true);
        if (compacted > 0) {
            stats_.memoryCompactions++;
        }
    }
}

double FragmentReservoir::getFragmentation() const {
    if (activePool_) {
        return activePool_->getFragmentation();
    }
    
    // Calculate our own fragmentation
    int totalSlots = instances_.size();
    int activeSlots = 0;
    for (const auto& inst : instances_) {
        if (inst.instanceId >= 0) {
            activeSlots++;
        }
    }
    
    if (totalSlots == 0) return 0.0;
    return 1.0 - static_cast<double>(activeSlots) / totalSlots;
}

bool FragmentReservoir::shouldCompact() const {
    return config_.autoCompact && getFragmentation() > config_.compactThreshold;
}

// ============================================================================
// Validation
// ============================================================================

bool FragmentReservoir::validate() const {
    // Check instance consistency
    for (size_t i = 0; i < instances_.size(); ++i) {
        const auto& inst = instances_[i];
        if (inst.instanceId < 0) continue;  // Free slot
        
        if (inst.isActive && inst.isGhost) {
            std::cerr << "Instance " << i << " is both active and ghost!" << std::endl;
            return false;
        }
        
        if (!isValidTemplateId(inst.templateId)) {
            std::cerr << "Instance " << i << " has invalid template ID!" << std::endl;
            return false;
        }
    }
    
    // Check index consistency
    for (const auto& [type, indices] : activeByType_) {
        for (int idx : indices) {
            if (!instances_[idx].isActive || instances_[idx].isGhost) {
                std::cerr << "Active index " << idx << " is not actually active!" << std::endl;
                return false;
            }
        }
    }
    
    for (int idx : ghostIndices_) {
        if (!instances_[idx].isGhost) {
            std::cerr << "Ghost index " << idx << " is not actually a ghost!" << std::endl;
            return false;
        }
    }
    
    return true;
}

// ============================================================================
// Missing Function Implementations
// ============================================================================

std::vector<int> FragmentReservoir::findInstancesInBox(const Vector3& min, const Vector3& max) const {
    std::vector<int> result;
    
    for (size_t i = 0; i < instances_.size(); ++i) {
        const auto& inst = instances_[i];
        if (!inst.isActive || inst.isGhost) continue;
        
        const Vector3& pos = inst.centerOfMass;
        if (pos.x >= min.x && pos.x <= max.x &&
            pos.y >= min.y && pos.y <= max.y &&
            pos.z >= min.z && pos.z <= max.z) {
            result.push_back(i);
        }
    }
    
    return result;
}

std::vector<int> FragmentReservoir::createMultipleInstances(int templateId, 
                                                           const std::vector<Vector3>& positions) {
    std::vector<int> result;
    result.reserve(positions.size());
    
    for (const auto& pos : positions) {
        int id = createInstance(templateId, pos);
        if (id >= 0) {
            result.push_back(id);
        }
    }
    
    return result;
}

int FragmentReservoir::deleteMultipleInstances(const std::vector<int>& instanceIds) {
    int deleted = 0;
    
    for (int id : instanceIds) {
        if (deleteInstance(id)) {
            deleted++;
        }
    }
    
    return deleted;
}

void FragmentReservoir::updateNeighborLists(int instanceId, 
                                           const std::vector<int>& proteinAtoms,
                                           const std::vector<int>& otherFragments) {
    if (!isValidInstanceId(instanceId)) return;
    
    FragmentInstance& inst = instances_[instanceId];
    inst.proteinNeighbors = proteinAtoms;
    inst.fragmentNeighbors = otherFragments;
    inst.lastNeighborUpdate = currentMCStep_;
}

const std::vector<int>& FragmentReservoir::getProteinNeighbors(int instanceId) const {
    static const std::vector<int> empty;
    
    if (!isValidInstanceId(instanceId)) {
        return empty;
    }
    
    return instances_[instanceId].proteinNeighbors;
}

const std::vector<int>& FragmentReservoir::getFragmentNeighbors(int instanceId) const {
    static const std::vector<int> empty;
    
    if (!isValidInstanceId(instanceId)) {
        return empty;
    }
    
    return instances_[instanceId].fragmentNeighbors;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc