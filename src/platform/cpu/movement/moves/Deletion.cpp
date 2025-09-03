#include "Deletion.hpp"
#include "../pool/ActivePool.hpp"
#include "../bias/CavityBias.hpp"
#include "../bias/CavityBiasCore.hpp"  // New cavity bias implementation
#include "../bias/UnifiedAcceptance.hpp"  // Unified acceptance probability
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../../../simulation/simulation.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

DeletionMove::DeletionMove(ActivePool* activePool, CavityManager* cavityManager, 
                           EnergyInterface* energyCalc, CavityBiasCore* cavityCore)
    : activePool_(activePool),
      cavityManager_(cavityManager),
      cavityCore_(cavityCore),
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
    
    // Calculate system volume from box dimensions
    MovementParams paramsWithVolume = params;
    if (paramsWithVolume.volumeNm3 <= 0.0) {
        paramsWithVolume.volumeNm3 = state.info.box[0] * state.info.box[1] * state.info.box[2];
    }
    
    // Calculate cavity bias for deletion if enabled
    // IMPORTANT: Calculate cavity bias in the post-deletion state (residue inactive)
    // This represents the reverse insertion probability
    double cavityBias = 1.0;
    double Vbox = paramsWithVolume.volumeNm3;
    double Vcav_after = Vbox;
    
    if (paramsWithVolume.useCavityBias && cavityCore_) {
        // Use new CavityBiasCore - calculate after-state cavity volume
        // Prefer FAST_APPROX here to reduce discretization bias in DB tests
        cavityCore_->invalidateCache();
        CavityMode mode = CavityMode::FAST_APPROX;
        Vcav_after = std::max(1e-30, cavityCore_->calculateCavityVolume(state, mode));
        cavityBias = Vcav_after / Vbox;
    } else if (paramsWithVolume.useCavityBias && cavityManager_) {
        // Fallback to legacy CavityManager
        cavityManager_->invalidateCache();
        cavityBias = cavityManager_->calculateCavityBiasFactor(state);
        Vcav_after = cavityBias * Vbox;
    }
    
    // Restore active flag before acceptance decision
    state.residues[targetResIdx].active = true;
    
    // Calculate deletion acceptance probability (match legacy LogSpace path)
    double acceptProb;
    bool useCavityBiasFlag = paramsWithVolume.useCavityBias && (cavityCore_ || cavityManager_);
    if (useCavityBiasFlag && paramsWithVolume.thermalLambdaNm != 1.0) {
        // With cavity + Lambda³
        acceptProb = utils::LogSpaceCalculator::calculateDeletionProbabilityWithCavityAndLambda(
            n_before,
            deltaE,
            paramsWithVolume.beta,
            paramsWithVolume.chemicalPotential,
            cavityBias,                // fraction Vcav_after/V
            paramsWithVolume.volumeNm3,
            paramsWithVolume.thermalLambdaNm,
            paramsWithVolume.useLogSpace
        );
    } else if (useCavityBiasFlag) {
        // With cavity only
        acceptProb = utils::LogSpaceCalculator::calculateDeletionProbabilityWithCavity(
            n_before,
            deltaE,
            paramsWithVolume.beta,
            paramsWithVolume.chemicalPotential,
            cavityBias,                // fraction Vcav_after/V
            paramsWithVolume.volumeNm3,
            paramsWithVolume.useLogSpace
        );
    } else {
        // No cavity
        acceptProb = calculateDeletionProbability(
            n_before,
            deltaE,
            paramsWithVolume
        );
    }
    result.acceptanceProbability = acceptProb;
    result.cavityBiasFactor = cavityBias;  // Store cavity bias factor in result
    
    // Accept or reject
    bool accepted = utils::RandomUtils::metropolisAccept(acceptProb);
    result.accepted = accepted;
    
    if (accepted) {
        // Permanently delete the residue
        // Mark as inactive again for deletion
        state.residues[targetResIdx].active = false;
        
        // Directly manipulate state to ensure correct deletion
        // This ensures we delete exactly the residue we tested
        state.removeResidue(targetResIdx);
        
        // Sync pool with the updated state
        if (activePool_) {
            activePool_->syncFromState(state);
        }
        
        // Invalidate cavity cache after accepted deletion
        if (cavityCore_) cavityCore_->invalidateCache();
        if (cavityManager_) cavityManager_->invalidateCache();
        
        stats_.acceptedDeletions++;
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
    
    // Use version with thermal wavelength if specified
    if (params.thermalLambdaNm != 1.0) {
        // Use cavity bias = 1.0 for non-cavity case
        return utils::LogSpaceCalculator::calculateDeletionProbabilityWithCavityAndLambda(
            n,
            deltaE,
            params.beta,
            params.chemicalPotential,
            1.0,  // cavity bias = 1.0 when not using cavity
            volumeNm3,
            params.thermalLambdaNm,
            params.useLogSpace
        );
    } else {
        return utils::LogSpaceCalculator::calculateDeletionProbability(
            n,
            deltaE,
            params.beta,
            params.chemicalPotential,
            volumeNm3,
            params.useLogSpace
        );
    }
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
