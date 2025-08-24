#include "Deletion.hpp"
#include "../pool/ActivePool.hpp"
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../../../simulation/simulation.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

DeletionMove::DeletionMove(ActivePool* activePool, EnergyInterface* energyCalc)
    : activePool_(activePool),
      energyCalc_(energyCalc) {
    resetStatistics();
}

DeletionMove::~DeletionMove() = default;

MovementResult DeletionMove::attemptDeletion(MCState& state, const MovementParams& params) {
    return performDeletion(state, params, -1);
}

MovementResult DeletionMove::performDeletion(MCState& state, const MovementParams& params, int residueIndex) {
    MovementResult result;
    result.moveType = "delete";
    
    // Check for valid box dimensions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        result.accepted = false;
        result.rejectReason = "Invalid box dimensions";
        result.energyChange = 0.0;
        result.acceptanceProbability = 0.0;
        stats_.totalAttempts++;
        return result;
    }
    
    // Check if there are any molecules to delete
    if (state.activeResidueCount == 0) {
        result.accepted = false;
        result.rejectReason = "No molecules to delete";
        stats_.totalAttempts++;
        stats_.rejectedEmpty++;
        return result;
    }
    
    // Select residue to delete
    int targetResIdx = residueIndex;
    if (targetResIdx < 0) {
        targetResIdx = selectResidueForDeletion(state);
    }
    
    if (targetResIdx < 0 || targetResIdx >= state.activeResidueCount) {
        result.accepted = false;
        result.rejectReason = "Invalid residue index";
        stats_.totalAttempts++;
        return result;
    }
    
    result.residueIndex = targetResIdx;
    
    // Calculate energy before deletion
    int n_before = state.activeResidueCount;  // Store n before deletion
    simulation::Simulation::computeSystemEnergyPBCCutoff(state);
    double energyBefore = 0.0;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        energyBefore += state.residues[i].energy_vdw;
        energyBefore += state.residues[i].energy_elec;
    }
    energyBefore *= 0.5;  // Account for double counting
    
    // Save residue information for potential restoration
    MCResidue savedResidue = state.residues[targetResIdx];
    std::vector<MCAtom> savedAtoms;
    for (int i = 0; i < savedResidue.atomCount; ++i) {
        int atomIdx = savedResidue.atomStart + i;
        if (atomIdx < state.activeAtomCount) {
            savedAtoms.push_back(state.atoms[atomIdx]);
        }
    }
    
    // Temporarily mark residue as inactive
    state.residues[targetResIdx].active = false;
    
    // Calculate energy after deletion (with residue marked inactive)
    simulation::Simulation::computeSystemEnergyPBCCutoff(state);
    double energyAfter = 0.0;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (i != targetResIdx && state.residues[i].active) {
            energyAfter += state.residues[i].energy_vdw;
            energyAfter += state.residues[i].energy_elec;
        }
    }
    energyAfter *= 0.5;
    
    double deltaE = energyAfter - energyBefore;
    result.energyChange = deltaE;
    
    // Restore active flag before acceptance decision
    state.residues[targetResIdx].active = true;
    
    // Calculate system volume from box dimensions
    MovementParams paramsWithVolume = params;
    if (paramsWithVolume.volumeNm3 <= 0.0) {
        paramsWithVolume.volumeNm3 = state.info.box[0] * state.info.box[1] * state.info.box[2];
    }
    
    // Calculate deletion acceptance probability (use n_before)
    double acceptProb = calculateDeletionProbability(
        n_before,  // n before deletion
        deltaE,
        paramsWithVolume
    );
    result.acceptanceProbability = acceptProb;
    
    // Accept or reject
    bool accepted = utils::RandomUtils::metropolisAccept(acceptProb);
    result.accepted = accepted;
    
    if (accepted) {
        // Permanently delete the residue
        // Mark as inactive again for deletion
        state.residues[targetResIdx].active = false;
        
        // First, update active pool
        auto activeIndices = activePool_->getActiveResidueIndices();
        if (targetResIdx < static_cast<int>(activeIndices.size())) {
            int poolResIdx = activeIndices[targetResIdx];
            if (activePool_->deleteResidue(poolResIdx)) {
                // Sync state with pool
                activePool_->syncToState(state);
                stats_.acceptedDeletions++;
            } else {
                // Pool deletion failed, restore
                state.residues[targetResIdx].active = true;
                result.accepted = false;
                result.rejectReason = "Active pool deletion failed";
            }
        } else {
            // Direct state manipulation (fallback)
            state.removeResidue(targetResIdx);
            stats_.acceptedDeletions++;
        }
    }
    // else: residue is already active from restoration before acceptance decision
    
    // Update statistics
    stats_.totalAttempts++;
    updateStatistics(accepted, deltaE);
    
    return result;
}

int DeletionMove::selectResidueForDeletion(const MCState& state) {
    if (state.activeResidueCount == 0) {
        return -1;
    }
    
    // Get list of active residues
    std::vector<int> activeResidues;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (state.residues[i].active) {
            activeResidues.push_back(i);
        }
    }
    
    if (activeResidues.empty()) {
        return -1;
    }
    
    // Select random active residue
    int idx = utils::RandomUtils::uniformInt(0, static_cast<int>(activeResidues.size()) - 1);
    return activeResidues[idx];
}

double DeletionMove::calculateDeletionProbability(
    int n,
    double deltaE,
    const MovementParams& params) {
    
    // Get system volume in nm^3
    double volumeNm3 = params.volumeNm3;
    if (volumeNm3 <= 0.0) {
        // Note: Caller should calculate volume from state.info.box
        // before calling this function
        volumeNm3 = 1.0;  // Default fallback
    }
    
    return utils::LogSpaceCalculator::calculateDeletionProbability(
        n,
        deltaE,
        params.beta,
        params.chemicalPotential,
        volumeNm3,
        params.useLogSpace
    );
}

double DeletionMove::calculateResidueEnergy(const MCState& state, int residueIndex) {
    if (residueIndex < 0 || residueIndex >= state.activeResidueCount) {
        return 0.0;
    }
    
    const MCResidue& residue = state.residues[residueIndex];
    if (!residue.active) {
        return 0.0;
    }
    
    // Return pre-calculated energy if available
    return residue.energy_vdw + residue.energy_elec;
}

void DeletionMove::resetStatistics() {
    stats_ = Statistics();
}

void DeletionMove::updateStatistics(bool accepted, double energyChange) {
    if (accepted) {
        stats_.averageEnergyChange = (stats_.averageEnergyChange * stats_.acceptedDeletions + energyChange) / 
                                     (stats_.acceptedDeletions + 1);
    }
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc