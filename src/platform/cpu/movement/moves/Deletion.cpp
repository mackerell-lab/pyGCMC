#include "Deletion.hpp"
#include "../../energy/EnergyModule.hpp"
#include "../pool/ActivePool.hpp"
#include "../bias/CavityBias.hpp"
#include "../bias/CavityBiasCore.hpp"  // New cavity bias implementation
#include "../bias/UnifiedAcceptance.hpp"  // Unified acceptance probability
#include "../common/MoveCommon.hpp"
#include "../common/MovementUtils.hpp"
#include "../gcmc/GCMCAcceptance.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include <limits>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;
using gcmc::GCMCAcceptance;

namespace {

struct DeletionCbmcResult {
    double logWeight = 0.0;
    double weight = 1.0;
    int trialsUsed = 1;
    bool computed = false;
};

DeletionCbmcResult computeDeletionCbmcWeight(
    MCState& state,
    int residueIndex,
    const MovementParams& params,
    const std::vector<MCAtom>& savedAtoms) {
    DeletionCbmcResult result;
    if (!params.useConfigBias || !params.useConfigBiasForInsertion || params.numConfigTrials <= 1) {
        return result;
    }
    const int numTrials = std::max(params.numConfigTrials, 1);
    if (savedAtoms.empty() || residueIndex < 0 || residueIndex >= state.activeResidueCount) {
        return result;
    }

    const MCResidue& residue = state.residues[residueIndex];
    const int atomCount = std::min(residue.atomCount, static_cast<int>(savedAtoms.size()));
    if (atomCount <= 0) {
        return result;
    }

    std::vector<Vector3> originalPositions;
    originalPositions.reserve(atomCount);
    for (int i = 0; i < atomCount; ++i) {
        originalPositions.emplace_back(savedAtoms[i].x, savedAtoms[i].y, savedAtoms[i].z);
    }

    Vector3 center(0.0, 0.0, 0.0);
    for (const auto& pos : originalPositions) {
        center.x += pos.x;
        center.y += pos.y;
        center.z += pos.z;
    }
    center.x /= static_cast<double>(atomCount);
    center.y /= static_cast<double>(atomCount);
    center.z /= static_cast<double>(atomCount);

    std::vector<Vector3> relativePositions = originalPositions;
    for (auto& pos : relativePositions) {
        pos = pos - center;
    }

    auto applyPositions = [&](const std::vector<Vector3>& positions) {
        for (int i = 0; i < atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < 0 || atomIdx >= state.activeAtomCount) continue;
            state.atoms[atomIdx].x = static_cast<float>(positions[i].x);
            state.atoms[atomIdx].y = static_cast<float>(positions[i].y);
            state.atoms[atomIdx].z = static_cast<float>(positions[i].z);
            if (i < static_cast<int>(state.residues[residueIndex].atoms.size())) {
                auto& atom = state.residues[residueIndex].atoms[i];
                atom.x = state.atoms[atomIdx].x;
                atom.y = state.atoms[atomIdx].y;
                atom.z = state.atoms[atomIdx].z;
            }
        }
    };

    auto restoreOriginal = [&]() {
        applyPositions(originalPositions);
    };

    auto computeResidueEnergy = [&]() -> double {
        platform::cpu::computeResidueEnergyCutoffPBC(state, residueIndex);
        const auto& res = state.residues[residueIndex];
        return static_cast<double>(res.energy_vdw + res.energy_elec);
    };

    restoreOriginal();

    std::vector<double> logWeights;
    logWeights.reserve(numTrials);

    double energy = computeResidueEnergy();
    if (std::isfinite(energy)) {
        logWeights.push_back(-params.beta * energy);
    }

    Vector3 box(state.info.box[0], state.info.box[1], state.info.box[2]);
    std::vector<Vector3> trialPositions(atomCount);

    for (int trial = 1; trial < numTrials; ++trial) {
        Quaternion quat = utils::RotationUtils::generateRandomQuaternion();
        double matrix[3][3];
        utils::RotationUtils::quaternionToMatrix(quat, matrix);
        Vector3 translation = center;
        if (params.configTranslationRange > 0.0) {
            translation = translation + utils::RandomUtils::randomVector(params.configTranslationRange);
        }

        for (int i = 0; i < atomCount; ++i) {
            Vector3 rotated = utils::RotationUtils::rotateVector(relativePositions[i], matrix);
            Vector3 placed = translation + rotated;
            if (box.x > 0.0 && box.y > 0.0 && box.z > 0.0) {
                placed = utils::PBCUtils::applyPBC(placed, box);
            }
            trialPositions[i] = placed;
        }

        applyPositions(trialPositions);
        double trialEnergy = computeResidueEnergy();
        if (std::isfinite(trialEnergy)) {
            logWeights.push_back(-params.beta * trialEnergy);
        }
    }

    restoreOriginal();

    const int keff = static_cast<int>(logWeights.size());
    if (keff == 0) {
        return result;
    }

    double logW = utils::LogSpaceCalculator::logSumExp(logWeights);
    double logK = std::log(static_cast<double>(std::max(keff, 1)));
    result.logWeight = logW - logK;
    result.weight = std::exp(result.logWeight);
    result.trialsUsed = keff;
    result.computed = true;
    return result;
}

void finalizeDeletionAcceptance(
    MovementResult& result,
    const MovementParams& params,
    int moleculeType,
    int countBefore) {
    result.hasGrandTerms = true;

    result.grandTerms.speciesId = moleculeType;
    result.grandTerms.countBefore = countBefore;
    result.grandTerms.countAfter = std::max(countBefore - 1, 0);
    result.grandTerms.beta = params.beta;
    result.grandTerms.chemicalPotential = params.chemicalPotential;
    result.grandTerms.deltaEnergy = result.energyChange;
    result.grandTerms.logVolume = result.logVolume;
    result.grandTerms.logLambda3 = result.logLambda3;
    result.grandTerms.logProposalForward = result.logProposalForward;
    result.grandTerms.logProposalReverse = result.logProposalReverse;
    result.grandTerms.logCavityForward = 0.0;
    result.grandTerms.logCavityReverse = result.logCavityFactor;
    result.grandTerms.logRosenbluthForward = 0.0;
    result.grandTerms.logRosenbluthReverse = result.logWReverse;
    result.grandTerms.logExtraForward = 0.0;
    result.grandTerms.logExtraReverse = 0.0;

    result.grandEvaluation = GCMCAcceptance::evaluate(
        result.grandTerms,
        GCMCAcceptance::MoveType::DELETION);
    result.logAcceptanceRatio = result.grandEvaluation.logRatio;
    result.acceptanceProbability = result.grandEvaluation.probability;
}

inline double getStoredRosenbluthWeight(const MCState& state, int residueIndex) {
    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        return state.residues[residueIndex].cbmcInsertionWeight;
    }
    return 1.0;
}

} // namespace

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

    result.lambdaNm = (params.thermalLambdaNm > 0.0) ? params.thermalLambdaNm : 1.0;
    result.logLambda3 = (std::abs(result.lambdaNm - 1.0) > 1e-12)
        ? 3.0 * std::log(result.lambdaNm)
        : 0.0;
    
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
    result.moleculeType = state.residues[targetResIdx].type;
    const int speciesCountBefore = move_common::countActiveResiduesOfType(state, result.moleculeType);
    const auto scheduler = params.getBiasedMoveProbabilitySet(result.moleculeType, speciesCountBefore);
    result.logProposalForward = utils::safeLogProbability(scheduler.deletion);
    result.logProposalReverse = utils::safeLogProbability(scheduler.insertion);
    
    // Calculate energy before deletion
    platform::cpu::computeSystemEnergyPBCCutoff(state);
    double energyBefore = move_common::sumResiduePairEnergy(state);
    
    // Save residue information for potential restoration
    MCResidue savedResidue = state.residues[targetResIdx];
    std::vector<MCAtom> savedAtoms;
    for (int i = 0; i < savedResidue.atomCount; ++i) {
        int atomIdx = savedResidue.atomStart + i;
        if (atomIdx < state.activeAtomCount) {
            savedAtoms.push_back(state.atoms[atomIdx]);
        }
    }

    double storedWeight = std::max(getStoredRosenbluthWeight(state, targetResIdx), 1e-30);
    DeletionCbmcResult cbmcInfo = computeDeletionCbmcWeight(state, targetResIdx, params, savedAtoms);
    if (cbmcInfo.computed) {
        result.cbmcTrialsUsed = cbmcInfo.trialsUsed;
        result.rosenbluthWeight = cbmcInfo.weight;
        result.logWReverse = cbmcInfo.logWeight;
    } else {
        result.cbmcTrialsUsed = 1;
        result.rosenbluthWeight = storedWeight;
        result.logWReverse = utils::safeLogProbability(storedWeight);
    }
    result.logWForward = 0.0;
    
    // Temporarily mark residue as inactive
    state.residues[targetResIdx].active = false;
    
    // Calculate energy after deletion (with residue marked inactive)
    platform::cpu::computeSystemEnergyPBCCutoff(state);
    double energyAfter = move_common::sumResiduePairEnergy(state, targetResIdx, true);
    
    double deltaE = energyAfter - energyBefore;
    result.energyChange = deltaE;
    
    // Calculate system volume from box dimensions
    MovementParams paramsWithVolume = params;
    if (paramsWithVolume.volumeNm3 <= 0.0) {
        paramsWithVolume.volumeNm3 = state.info.box[0] * state.info.box[1] * state.info.box[2];
    }
    result.volumeNm3 = paramsWithVolume.volumeNm3;
    result.logVolume = std::log(std::max(result.volumeNm3, 1e-30));
    
    // Calculate cavity bias for deletion if enabled
    // IMPORTANT: Calculate cavity bias in the post-deletion state (residue inactive)
    // This represents the reverse insertion probability
    double cavityBias = 1.0;
    double Vbox = paramsWithVolume.volumeNm3;
    double Vcav_after = Vbox;
    
    if (paramsWithVolume.useCavityBias && cavityCore_) {
        cavityCore_->invalidateCache();
        CavityMode mode = CavityMode::FAST_APPROX;
        Vcav_after = std::max(1e-30, cavityCore_->calculateCavityVolume(state, mode, result.moleculeType));
        cavityBias = Vcav_after / Vbox;
    } else if (paramsWithVolume.useCavityBias && cavityManager_) {
        cavityManager_->invalidateCache();
        cavityBias = cavityManager_->calculateCavityBiasFactor(state, result.moleculeType);
        Vcav_after = cavityBias * Vbox;
    }
    
    // Restore active flag before acceptance decision
    state.residues[targetResIdx].active = true;

    result.cavityBiasFactor = cavityBias;
    bool useCavityBiasFlag = paramsWithVolume.useCavityBias && (cavityCore_ || cavityManager_);
    const double populationNorm = std::max(1, speciesCountBefore);
    if (useCavityBiasFlag) {
        result.logCavityFactor = utils::safeLogProbability(cavityBias);
        result.cavityVolumeNm3 = Vcav_after;
        result.effectiveVolumeNm3 = Vcav_after / populationNorm;
    } else {
        result.logCavityFactor = 0.0;
        result.cavityVolumeNm3 = result.volumeNm3;
        result.effectiveVolumeNm3 = result.volumeNm3 / populationNorm;
    }

    finalizeDeletionAcceptance(result, paramsWithVolume, result.moleculeType, speciesCountBefore);
    
    // Accept or reject
    bool accepted = utils::RandomUtils::metropolisAccept(result.acceptanceProbability);
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
