#include "ActivePool.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include <algorithm>
#include <numeric>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

ActivePool::ActivePool(int maxAtoms, int maxResidues)
    : maxAtoms_(maxAtoms),
      maxResidues_(maxResidues),
      fragmentationThreshold_(0.3),
      activeAtomCount_(0),
      activeResidueCount_(0),
      nextFreeAtomIndex_(0),
      nextFreeResidueIndex_(0) {
    
    // Pre-allocate memory
    atoms_.reserve(maxAtoms);
    residueInfo_.resize(maxResidues);
    atomActive_.resize(maxAtoms, false);
    residueActive_.resize(maxResidues, false);
    
    // Initialize free residue slots
    for (int i = 0; i < maxResidues; ++i) {
        freeResidueSlots_.insert(i);
    }
    
    // Reset statistics
    resetStatistics();
}

ActivePool::~ActivePool() = default;

int ActivePool::insertMolecule(const std::vector<MCAtom>& atoms, int resType) {
    // Check if we have space
    if (!canInsert(static_cast<int>(atoms.size()))) {
        return -1;
    }
    
    // Find a free residue slot
    int resIdx = findFreeResidueSlot();
    if (resIdx < 0) {
        return -1;
    }
    
    // Find space for atoms (append-only strategy for now)
    int atomStartIdx = nextFreeAtomIndex_;
    if (atomStartIdx + static_cast<int>(atoms.size()) > maxAtoms_) {
        // Try to compact if fragmented
        if (shouldCompact()) {
            compact(true);
            atomStartIdx = nextFreeAtomIndex_;
            if (atomStartIdx + static_cast<int>(atoms.size()) > maxAtoms_) {
                return -1;  // Still no space after compaction
            }
        } else {
            return -1;
        }
    }
    
    // Copy atoms to pool
    copyAtomsToPool(atoms, atomStartIdx);
    
    // Update residue info
    ResidueMetadata& resInfo = residueInfo_[resIdx];
    resInfo.residueType = resType;
    resInfo.active = true;
    resInfo.atomStartIndex = atomStartIdx;
    resInfo.atomCount = static_cast<int>(atoms.size());
    resInfo.insertionTime = static_cast<double>(stats_.totalInserts);  // Use insert count as pseudo-time
    
    // Update active masks
    residueActive_[resIdx] = true;
    for (int i = 0; i < resInfo.atomCount; ++i) {
        atomActive_[atomStartIdx + i] = true;
    }
    
    // Update counters
    activeAtomCount_ += resInfo.atomCount;
    activeResidueCount_++;
    nextFreeAtomIndex_ += resInfo.atomCount;
    freeResidueSlots_.erase(resIdx);
    
    // Update statistics
    stats_.totalInserts++;
    stats_.peakAtoms = std::max(stats_.peakAtoms, activeAtomCount_);
    stats_.peakResidues = std::max(stats_.peakResidues, activeResidueCount_);
    
    // Update center of mass
    updateCenterOfMass(resIdx);
    
    return resIdx;
}

bool ActivePool::deleteResidue(int resIdx) {
    if (resIdx < 0 || resIdx >= maxResidues_ || !residueActive_[resIdx]) {
        return false;
    }
    
    // Mark residue as inactive
    residueActive_[resIdx] = false;
    ResidueMetadata& resInfo = residueInfo_[resIdx];
    resInfo.active = false;
    
    // Mark atoms as inactive
    for (int i = 0; i < resInfo.atomCount; ++i) {
        atomActive_[resInfo.atomStartIndex + i] = false;
    }
    
    // Update counters
    activeAtomCount_ -= resInfo.atomCount;
    activeResidueCount_--;
    freeResidueSlots_.insert(resIdx);
    
    // Update statistics
    stats_.totalDeletes++;
    updateFragmentationStats();
    
    return true;
}

int ActivePool::compact(bool force) {
    if (!force && !shouldCompact()) {
        return 0;
    }
    
    int compactedAtoms = 0;
    int writeIdx = 0;
    
    // Create new atom array with only active atoms
    std::vector<MCAtom> compactedAtoms_;
    compactedAtoms_.reserve(activeAtomCount_);
    
    // Compact atoms and update residue indices
    for (int resIdx = 0; resIdx < maxResidues_; ++resIdx) {
        if (!residueActive_[resIdx]) continue;
        
        ResidueMetadata& resInfo = residueInfo_[resIdx];
        int oldStart = resInfo.atomStartIndex;
        int newStart = writeIdx;
        
        // Copy active atoms
        for (int i = 0; i < resInfo.atomCount; ++i) {
            if (atomActive_[oldStart + i]) {
                compactedAtoms_.push_back(atoms_[oldStart + i]);
                writeIdx++;
            }
        }
        
        // Update residue atom indices
        resInfo.atomStartIndex = newStart;
        compactedAtoms += (oldStart - newStart);
    }
    
    // Replace atom array
    atoms_ = std::move(compactedAtoms_);
    nextFreeAtomIndex_ = writeIdx;
    
    // Reset active masks for compacted atoms
    std::fill(atomActive_.begin(), atomActive_.end(), false);
    for (int i = 0; i < writeIdx; ++i) {
        atomActive_[i] = true;
    }
    
    // Update statistics
    stats_.compactions++;
    updateFragmentationStats();
    
    return compactedAtoms;
}

void ActivePool::queueInsert(const std::vector<MCAtom>& atoms, int resType) {
    insertQueue_.push({atoms, resType});
}

void ActivePool::queueDelete(int resIdx) {
    deleteQueue_.push(resIdx);
}

std::pair<int, int> ActivePool::flushBatch() {
    int inserted = 0;
    int deleted = 0;
    
    // Process deletions first (to free space)
    while (!deleteQueue_.empty()) {
        if (deleteResidue(deleteQueue_.front())) {
            deleted++;
        }
        deleteQueue_.pop();
    }
    
    // Process insertions
    while (!insertQueue_.empty()) {
        const auto& op = insertQueue_.front();
        if (insertMolecule(op.atoms, op.resType) >= 0) {
            inserted++;
        }
        insertQueue_.pop();
    }
    
    stats_.batchOperations++;
    
    return {inserted, deleted};
}

void ActivePool::syncToState(MCState& state) {
    // Clear current state atoms and residues
    state.atoms.clear();
    state.residues.clear();
    state.atoms.reserve(activeAtomCount_);
    state.residues.reserve(activeResidueCount_);
    
    // Copy active atoms and build residues
    for (int resIdx = 0; resIdx < maxResidues_; ++resIdx) {
        if (!residueActive_[resIdx]) continue;
        
        const ResidueMetadata& resInfo = residueInfo_[resIdx];
        
        // Create MCResidue
        MCResidue mcRes;
        mcRes.atomStart = static_cast<int>(state.atoms.size());
        mcRes.atomCount = resInfo.atomCount;
        mcRes.type = resInfo.residueType;
        mcRes.active = true;
        
        // Copy atoms for this residue
        for (int i = 0; i < resInfo.atomCount; ++i) {
            int atomIdx = resInfo.atomStartIndex + i;
            if (atomIdx < static_cast<int>(atoms_.size())) {
                state.atoms.push_back(atoms_[atomIdx]);
            }
        }
        
        state.residues.push_back(mcRes);
    }
    
    // Update state counters
    state.activeAtomCount = static_cast<int>(state.atoms.size());
    state.activeResidueCount = static_cast<int>(state.residues.size());
}

void ActivePool::syncFromState(const MCState& state) {
    // Clear current pool
    atoms_.clear();
    std::fill(atomActive_.begin(), atomActive_.end(), false);
    std::fill(residueActive_.begin(), residueActive_.end(), false);
    freeResidueSlots_.clear();
    
    // Reset counters
    activeAtomCount_ = 0;
    activeResidueCount_ = 0;
    nextFreeAtomIndex_ = 0;
    nextFreeResidueIndex_ = 0;
    
    // Copy atoms from state
    atoms_ = state.atoms;
    activeAtomCount_ = state.activeAtomCount;
    nextFreeAtomIndex_ = activeAtomCount_;
    
    // Rebuild residue info
    for (int i = 0; i < state.activeResidueCount; ++i) {
        const MCResidue& mcRes = state.residues[i];
        if (!mcRes.active) continue;
        
        ResidueMetadata& resInfo = residueInfo_[i];
        resInfo.residueType = mcRes.type;
        resInfo.active = true;
        resInfo.atomStartIndex = mcRes.atomStart;
        resInfo.atomCount = mcRes.atomCount;
        
        residueActive_[i] = true;
        for (int j = 0; j < mcRes.atomCount; ++j) {
            atomActive_[mcRes.atomStart + j] = true;
        }
        
        updateCenterOfMass(i);
        activeResidueCount_++;
    }
    
    // Update free slots
    for (int i = activeResidueCount_; i < maxResidues_; ++i) {
        freeResidueSlots_.insert(i);
    }
    nextFreeResidueIndex_ = activeResidueCount_;
}

double ActivePool::getFragmentation() const {
    if (nextFreeAtomIndex_ == 0) return 0.0;
    return 1.0 - static_cast<double>(activeAtomCount_) / nextFreeAtomIndex_;
}

std::pair<int, int> ActivePool::getActiveCounts() const {
    return {activeAtomCount_, activeResidueCount_};
}

std::vector<int> ActivePool::getActiveResidueIndices() const {
    std::vector<int> indices;
    indices.reserve(activeResidueCount_);
    
    for (int i = 0; i < maxResidues_; ++i) {
        if (residueActive_[i]) {
            indices.push_back(i);
        }
    }
    
    return indices;
}

bool ActivePool::isResidueActive(int resIdx) const {
    return resIdx >= 0 && resIdx < maxResidues_ && residueActive_[resIdx];
}

ActivePool::ResidueMetadata* ActivePool::getResidueMetadata(int resIdx) {
    if (resIdx >= 0 && resIdx < maxResidues_) {
        return &residueInfo_[resIdx];
    }
    return nullptr;
}

const ActivePool::ResidueMetadata* ActivePool::getResidueMetadata(int resIdx) const {
    if (resIdx >= 0 && resIdx < maxResidues_) {
        return &residueInfo_[resIdx];
    }
    return nullptr;
}

void ActivePool::resetStatistics() {
    stats_ = Statistics();
}

bool ActivePool::canInsert(int atomCount) const {
    bool hasResidueSlot = !freeResidueSlots_.empty() || nextFreeResidueIndex_ < maxResidues_;
    bool hasAtomSpace = (nextFreeAtomIndex_ + atomCount <= maxAtoms_) || 
                        (shouldCompact() && activeAtomCount_ + atomCount <= maxAtoms_);
    return hasResidueSlot && hasAtomSpace;
}

// Private helper functions

int ActivePool::findFreeResidueSlot() {
    if (!freeResidueSlots_.empty()) {
        return *freeResidueSlots_.begin();
    }
    
    if (nextFreeResidueIndex_ < maxResidues_) {
        return nextFreeResidueIndex_++;
    }
    
    return -1;
}

int ActivePool::findFreeAtomRange(int count) {
    // Simple append-only strategy for now
    if (nextFreeAtomIndex_ + count <= maxAtoms_) {
        return nextFreeAtomIndex_;
    }
    
    // Could implement more sophisticated free range tracking later
    return -1;
}

void ActivePool::updateFragmentationStats() {
    double frag = getFragmentation();
    stats_.averageFragmentation = (stats_.averageFragmentation * (stats_.totalInserts + stats_.totalDeletes - 1) + frag) / 
                                  (stats_.totalInserts + stats_.totalDeletes);
}

bool ActivePool::shouldCompact() const {
    return getFragmentation() > fragmentationThreshold_;
}

void ActivePool::copyAtomsToPool(const std::vector<MCAtom>& atoms, int startIdx) {
    // Ensure atoms_ vector has enough size
    while (static_cast<int>(atoms_.size()) < startIdx + static_cast<int>(atoms.size())) {
        atoms_.emplace_back();
    }
    
    // Copy atoms
    for (size_t i = 0; i < atoms.size(); ++i) {
        atoms_[startIdx + i] = atoms[i];
    }
}

void ActivePool::updateCenterOfMass(int resIdx) {
    if (!residueActive_[resIdx]) return;
    
    ResidueMetadata& resInfo = residueInfo_[resIdx];
    Vector3 com(0.0, 0.0, 0.0);
    
    for (int i = 0; i < resInfo.atomCount; ++i) {
        int atomIdx = resInfo.atomStartIndex + i;
        if (atomIdx < static_cast<int>(atoms_.size())) {
            const MCAtom& atom = atoms_[atomIdx];
            com.x += atom.x;
            com.y += atom.y;
            com.z += atom.z;
        }
    }
    
    if (resInfo.atomCount > 0) {
        com.x /= resInfo.atomCount;
        com.y /= resInfo.atomCount;
        com.z /= resInfo.atomCount;
    }
    
    resInfo.centerOfMass = com;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc