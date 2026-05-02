#include "Rotation.hpp"
#include "../../energy/EnergyModule.hpp"
#include "../pool/ActivePool.hpp"
#include "../bias/ConfigBias.hpp"
#include "../common/MoveCommon.hpp"
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

RotationMove::RotationMove(ActivePool* activePool,
                           ConfigBiasManager* configBiasManager,
                           EnergyInterface* energyCalc)
    : activePool_(activePool),
      configBiasManager_(configBiasManager),
      energyCalc_(energyCalc) {
    resetStatistics();
}

RotationMove::~RotationMove() = default;

MovementResult RotationMove::attemptRotation(MCState& state, const MovementParams& params) {
    if (params.useConfigBias && configBiasManager_) {
        return performConfigBiasRotation(state, params, -1);
    } else {
        return performSimpleRotation(state, params, -1);
    }
}

MovementResult RotationMove::performSimpleRotation(MCState& state, const MovementParams& params, int residueIndex) {
    MovementResult result;
    result.moveType = "rotate";

    // Check for valid box dimensions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        result.accepted = false;
        result.rejectReason = "Invalid box dimensions";
        result.energyChange = 0.0;
        result.acceptanceProbability = 0.0;
        stats_.totalAttempts++;
        return result;
    }

    // Check if there are any molecules to rotate
    if (state.activeResidueCount == 0) {
        result.accepted = false;
        result.rejectReason = "No molecules to rotate";
        stats_.totalAttempts++;
        stats_.rejectedEmpty++;
        return result;
    }

    // Select residue to rotate
    int targetResIdx = residueIndex;
    if (targetResIdx < 0) {
        targetResIdx = selectResidueForRotation(state);
    }

    if (targetResIdx < 0 || targetResIdx >= state.activeResidueCount) {
        result.accepted = false;
        result.rejectReason = "Invalid residue index";
        stats_.totalAttempts++;
        return result;
    }

    result.residueIndex = targetResIdx;
    int fragmentType = state.residues[targetResIdx].type;
    const auto scheduler = params.getMoveProbabilitySet(fragmentType);
    double logProb = utils::safeLogProbability(scheduler.rotation);
    result.logProposalForward = logProb;
    result.logProposalReverse = logProb;

    // Save original configuration
    RotationConfig originalConfig = saveConfiguration(state, targetResIdx);

    // Generate random rotation
    Quaternion rotation = generateRandomRotation(params.maxRotation);

    // Calculate energy before rotation
    platform::cpu::computeSystemEnergyPBCCutoff(state);
    double energyBefore = move_common::sumResiduePairEnergy(state);

    // Apply rotation around center of mass
    rotateResidue(state, targetResIdx, rotation);

    // Calculate energy after rotation
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
        // Restore original configuration
        restoreConfiguration(state, targetResIdx, originalConfig);
    } else {
        stats_.acceptedRotations++;
    }

    // Update statistics
    stats_.totalAttempts++;
    stats_.simpleRotations++;
    updateStatistics(accepted, false, params.maxRotation, deltaE, 1.0);

    return result;
}

MovementResult RotationMove::performConfigBiasRotation(MCState& state, const MovementParams& params, int residueIndex) {
    MovementResult result;
    result.moveType = "rotate";

    // Check if there are any molecules to rotate
    if (state.activeResidueCount == 0) {
        result.accepted = false;
        result.rejectReason = "No molecules to rotate";
        stats_.totalAttempts++;
        stats_.rejectedEmpty++;
        return result;
    }

    // Select residue to rotate
    int targetResIdx = residueIndex;
    if (targetResIdx < 0) {
        targetResIdx = selectResidueForRotation(state);
    }

    if (targetResIdx < 0 || targetResIdx >= state.activeResidueCount) {
        result.accepted = false;
        result.rejectReason = "Invalid residue index";
        stats_.totalAttempts++;
        return result;
    }

    result.residueIndex = targetResIdx;
    int fragmentType = state.residues[targetResIdx].type;
    const auto scheduler = params.getMoveProbabilitySet(fragmentType);
    double logProb = utils::safeLogProbability(scheduler.rotation);
    result.logProposalForward = logProb;
    result.logProposalReverse = logProb;

    // Perform configurational bias rotation
    ConfigBiasRotationResult biasResult = performConfigBiasRotationInternal(state, targetResIdx, params);

    result.accepted = biasResult.accepted;
    result.energyChange = biasResult.energyChange;
    result.configBiasFactor = biasResult.biasFactor;
    result.numConfigTrials = biasResult.totalConfigs;

    // Calculate final acceptance probability
    double acceptProb = std::min(1.0, std::exp(-params.beta * biasResult.energyChange) / biasResult.biasFactor);
    result.acceptanceProbability = acceptProb;

    // Update statistics
    stats_.totalAttempts++;
    stats_.configBiasRotations++;
    if (biasResult.accepted) {
        stats_.acceptedRotations++;
    }
    updateStatistics(biasResult.accepted, true, params.maxRotation, biasResult.energyChange, biasResult.biasFactor);

    return result;
}

RotationMove::ConfigBiasRotationResult RotationMove::performConfigBiasRotationInternal(
    MCState& state,
    int residueIndex,
    const MovementParams& params) {

    ConfigBiasRotationResult result;
    result.accepted = false;

    // Save original configuration
    RotationConfig originalConfig = saveConfiguration(state, residueIndex);

    // Generate and evaluate trial configurations
    std::vector<Configuration> configs;
    configs.reserve(params.numConfigTrials);

    for (int i = 0; i < params.numConfigTrials; ++i) {
        Configuration config;
        config.rotation = generateRandomRotation(params.maxRotation);
        config.index = i;

        // Apply rotation
        rotateResidue(state, residueIndex, config.rotation);

        // Calculate energy
        platform::cpu::computeSystemEnergyPBCCutoff(state);
        config.energy = move_common::sumResiduePairEnergy(state);

        configs.push_back(config);

        // Restore for next trial
        restoreConfiguration(state, residueIndex, originalConfig);
    }

    // Calculate Boltzmann probabilities
    configBiasManager_->calculateProbabilities(configs, params.beta, params.useLogSpace);

    // Select configuration based on probability
    int selectedIdx = configBiasManager_->selectByProbability(configs);

    // Apply selected configuration
    if (selectedIdx >= 0 && selectedIdx < static_cast<int>(configs.size())) {
        rotateResidue(state, residueIndex, configs[selectedIdx].rotation);

        // Calculate bias factor
        result.biasFactor = configBiasManager_->calculateBiasFactor(configs, selectedIdx);
        result.energyChange = configs[selectedIdx].energy - originalConfig.originalEnergy;
        result.selectedConfigIndex = selectedIdx;
        result.totalConfigs = params.numConfigTrials;

        // Accept with bias-corrected probability
        double acceptProb = std::min(1.0, std::exp(-params.beta * result.energyChange) / result.biasFactor);
        result.accepted = utils::RandomUtils::metropolisAccept(acceptProb);

        if (!result.accepted) {
            // Restore original if rejected
            restoreConfiguration(state, residueIndex, originalConfig);
        }
    }

    return result;
}

int RotationMove::selectResidueForRotation(const MCState& state) {
    return move_common::selectRandomActiveResidue(state);
}

Quaternion RotationMove::generateRandomRotation(double maxAngle) {
    // Generate random rotation with angle constraint
    // Generate random angle up to maxAngle
    double angle = utils::RandomUtils::uniform(0.0, maxAngle);

    // Generate random axis (unit vector)
    double theta = utils::RandomUtils::uniform(0.0, 2.0 * M_PI);
    double phi = std::acos(utils::RandomUtils::uniform(-1.0, 1.0));

    Vector3 axis(
        std::sin(phi) * std::cos(theta),
        std::sin(phi) * std::sin(theta),
        std::cos(phi)
    );

    // Create quaternion from angle-axis
    return utils::RotationUtils::quaternionFromAxisAngle(axis, angle);
}

RotationMove::RotationConfig RotationMove::saveConfiguration(const MCState& state, int residueIndex) {
    RotationConfig config;

    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        const MCResidue& residue = state.residues[residueIndex];

        // Save atom positions
        config.originalPositions.reserve(residue.atomCount);
        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                const MCAtom& atom = state.atoms[atomIdx];
                config.originalPositions.emplace_back(atom.x, atom.y, atom.z);
            }
        }

        // Calculate center of mass
        config.center = calculateCenterOfMass(state, residueIndex);

        // Calculate original energy
        platform::cpu::computeSystemEnergyPBCCutoff(const_cast<MCState&>(state));
        config.originalEnergy = move_common::sumResiduePairEnergy(state);
    }

    return config;
}

void RotationMove::restoreConfiguration(MCState& state, int residueIndex, const RotationConfig& config) {
    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        MCResidue& residue = state.residues[residueIndex];

        for (int i = 0; i < residue.atomCount && i < static_cast<int>(config.originalPositions.size()); ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                MCAtom& atom = state.atoms[atomIdx];
                atom.x = config.originalPositions[i].x;
                atom.y = config.originalPositions[i].y;
                atom.z = config.originalPositions[i].z;
            }
        }
    }
}

void RotationMove::rotateResidue(MCState& state, int residueIndex, const Quaternion& quaternion) {
    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        MCResidue& residue = state.residues[residueIndex];

        // Calculate center of mass
        Vector3 center = calculateCenterOfMass(state, residueIndex);

        // Convert quaternion to rotation matrix
        double matrix[3][3];
        utils::RotationUtils::quaternionToMatrix(quaternion, matrix);

        // Rotate each atom around center
        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                MCAtom& atom = state.atoms[atomIdx];

                // Translate to origin
                Vector3 pos(atom.x - center.x, atom.y - center.y, atom.z - center.z);

                // Apply rotation
                Vector3 rotated = utils::RotationUtils::rotateVector(pos, matrix);

                // Translate back and apply PBC
                atom.x = rotated.x + center.x;
                atom.y = rotated.y + center.y;
                atom.z = rotated.z + center.z;

                // Apply periodic boundary conditions
                utils::PBCUtils::applyPBC(atom.x, atom.y, atom.z, state.info.box);
            }
        }
    }
}

Vector3 RotationMove::calculateCenterOfMass(const MCState& state, int residueIndex) {
    Vector3 com(0.0, 0.0, 0.0);

    if (residueIndex >= 0 && residueIndex < state.activeResidueCount) {
        const MCResidue& residue = state.residues[residueIndex];
        int count = 0;

        for (int i = 0; i < residue.atomCount; ++i) {
            int atomIdx = residue.atomStart + i;
            if (atomIdx < state.activeAtomCount) {
                const MCAtom& atom = state.atoms[atomIdx];
                com.x += atom.x;
                com.y += atom.y;
                com.z += atom.z;
                count++;
            }
        }

        if (count > 0) {
            com.x /= count;
            com.y /= count;
            com.z /= count;
        }
    }

    return com;
}

void RotationMove::resetStatistics() {
    stats_ = Statistics();
}

void RotationMove::updateStatistics(bool accepted, bool usedConfigBias,
                                    double angle, double energyChange, double biasFactor) {
    if (accepted) {
        stats_.averageAngle = (stats_.averageAngle * stats_.acceptedRotations + angle) /
                             (stats_.acceptedRotations + 1);
        stats_.averageEnergyChange = (stats_.averageEnergyChange * stats_.acceptedRotations + energyChange) /
                                     (stats_.acceptedRotations + 1);
    }

    if (usedConfigBias) {
        stats_.averageBiasFactor = (stats_.averageBiasFactor * (stats_.configBiasRotations - 1) + biasFactor) /
                                  stats_.configBiasRotations;
    }
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
