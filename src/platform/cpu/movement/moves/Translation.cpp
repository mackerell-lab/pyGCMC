#include "Translation.hpp"
#include "../../energy/EnergyModule.hpp"
#include "../pool/ActivePool.hpp"
#include "../common/MoveCommon.hpp"
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

TranslationMove::TranslationMove(ActivePool* activePool, EnergyInterface* energyCalc)
    : activePool_(activePool),
      energyCalc_(energyCalc) {
    resetStatistics();
}

TranslationMove::~TranslationMove() = default;

MovementResult TranslationMove::attemptTranslation(MCState& state, const MovementParams& params) {
    return performTranslation(state, params, -1);
}

MovementResult TranslationMove::performTranslation(MCState& state, const MovementParams& params, int residueIndex) {
    MovementResult result;
    result.moveType = "translate";
    
    // Check for valid box dimensions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        result.accepted = false;
        result.rejectReason = "Invalid box dimensions";
        result.energyChange = 0.0;
        result.acceptanceProbability = 0.0;
        stats_.totalAttempts++;
        return result;
    }
    
    // Check if there are any molecules to translate
    if (state.activeResidueCount == 0) {
        result.accepted = false;
        result.rejectReason = "No molecules to translate";
        stats_.totalAttempts++;
        stats_.rejectedEmpty++;
        return result;
    }
    
    // Select residue to translate
    int targetResIdx = residueIndex;
    if (targetResIdx < 0) {
        targetResIdx = selectResidueForTranslation(state);
    }
    
    if (targetResIdx < 0 || targetResIdx >= state.activeResidueCount) {
        result.accepted = false;
        result.rejectReason = "Invalid residue index";
        stats_.totalAttempts++;
        return result;
    }
    
    result.residueIndex = targetResIdx;
    MCResidue& residue = state.residues[targetResIdx];
    const auto scheduler = params.getMoveProbabilitySet(residue.type);
    double logProb = utils::safeLogProbability(scheduler.translation);
    result.logProposalForward = logProb;
    result.logProposalReverse = logProb;
    
    // Save original atom positions
    std::vector<Vector3> originalPositions = saveAtomPositions(state, targetResIdx);
    
    // Generate random displacement
    Vector3 displacement = generateDisplacement(params.maxTranslation);
    
    // Calculate energy before translation
    platform::cpu::computeSystemEnergyPBCCutoff(state);
    double energyBefore = move_common::sumResiduePairEnergy(state);
    
    // Apply translation
    translateResidue(state, targetResIdx, displacement);
    
    // Apply PBC to ensure atoms stay in box
    for (int i = 0; i < residue.atomCount; ++i) {
        int atomIdx = residue.atomStart + i;
        if (atomIdx < state.activeAtomCount) {
            MCAtom& atom = state.atoms[atomIdx];
            Vector3 pos(atom.x, atom.y, atom.z);
            Vector3 box(state.info.box[0], state.info.box[1], state.info.box[2]);
            pos = utils::PBCUtils::applyPBC(pos, box);
            atom.x = pos.x;
            atom.y = pos.y;
            atom.z = pos.z;
        }
    }
    
    // Calculate energy after translation
    platform::cpu::computeSystemEnergyPBCCutoff(state);
    double energyAfter = move_common::sumResiduePairEnergy(state);
    
    double deltaE = energyAfter - energyBefore;
    result.energyChange = deltaE;
    
    // Calculate acceptance probability (Metropolis criterion)
    double acceptProb = std::min(1.0, std::exp(-params.beta * deltaE));
    result.acceptanceProbability = acceptProb;
    
    // Accept or reject
    bool accepted = utils::RandomUtils::metropolisAccept(acceptProb);
    result.accepted = accepted;
    
    if (!accepted) {
        // Restore original positions
        restoreAtomPositions(state, targetResIdx, originalPositions);
    } else {
        stats_.acceptedTranslations++;
    }
    
    // Update statistics
    stats_.totalAttempts++;
    updateStatistics(accepted, displacement.norm(), deltaE);
    
    return result;
}

std::vector<MovementResult> TranslationMove::performBatchTranslations(
    MCState& state, 
    const MovementParams& params, 
    int numAttempts) {
    
    std::vector<MovementResult> results;
    results.reserve(numAttempts);
    
    for (int i = 0; i < numAttempts; ++i) {
        results.push_back(performTranslation(state, params, -1));
    }
    
    return results;
}

int TranslationMove::selectResidueForTranslation(const MCState& state) {
    return move_common::selectRandomActiveResidue(state);
}

Vector3 TranslationMove::generateDisplacement(double maxTranslation) {
    // Generate random displacement in each direction
    return Vector3(
        utils::RandomUtils::uniform(-maxTranslation, maxTranslation),
        utils::RandomUtils::uniform(-maxTranslation, maxTranslation),
        utils::RandomUtils::uniform(-maxTranslation, maxTranslation)
    );
}

std::vector<Vector3> TranslationMove::saveAtomPositions(const MCState& state, int residueIndex) {
    std::vector<Vector3> positions;
    
    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        const MCResidue& residue = state.residues[residueIndex];
        positions.reserve(residue.atomCount);
        
        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                const MCAtom& atom = state.atoms[atomIdx];
                positions.emplace_back(atom.x, atom.y, atom.z);
            }
        }
    }
    
    return positions;
}

void TranslationMove::restoreAtomPositions(MCState& state, int residueIndex, const std::vector<Vector3>& positions) {
    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        MCResidue& residue = state.residues[residueIndex];
        
        for (int i = 0; i < residue.atomCount && i < static_cast<int>(positions.size()); ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                MCAtom& atom = state.atoms[atomIdx];
                atom.x = positions[i].x;
                atom.y = positions[i].y;
                atom.z = positions[i].z;
            }
        }
    }
}

void TranslationMove::translateResidue(MCState& state, int residueIndex, const Vector3& displacement) {
    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        MCResidue& residue = state.residues[residueIndex];
        
        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                MCAtom& atom = state.atoms[atomIdx];
                atom.x += displacement.x;
                atom.y += displacement.y;
                atom.z += displacement.z;
            }
        }
    }
}

std::pair<double, double> TranslationMove::calculateEnergyChange(
    MCState& state, 
    int residueIndex, 
    const Vector3& displacement) {
    
    // Save current positions
    auto originalPos = saveAtomPositions(state, residueIndex);
    
    // Calculate energy before
    platform::cpu::computeSystemEnergyPBCCutoff(state);
    double energyBefore = move_common::sumResiduePairEnergy(state);
    
    // Apply translation
    translateResidue(state, residueIndex, displacement);
    
    // Calculate energy after
    platform::cpu::computeSystemEnergyPBCCutoff(state);
    double energyAfter = move_common::sumResiduePairEnergy(state);
    
    // Restore positions
    restoreAtomPositions(state, residueIndex, originalPos);
    
    return {energyBefore, energyAfter};
}

void TranslationMove::resetStatistics() {
    stats_ = Statistics();
}

void TranslationMove::updateStatistics(bool accepted, double displacement, double energyChange) {
    if (accepted) {
        stats_.averageDisplacement = (stats_.averageDisplacement * stats_.acceptedTranslations + displacement) / 
                                     (stats_.acceptedTranslations + 1);
        stats_.averageEnergyChange = (stats_.averageEnergyChange * stats_.acceptedTranslations + energyChange) / 
                                     (stats_.acceptedTranslations + 1);
    }
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
