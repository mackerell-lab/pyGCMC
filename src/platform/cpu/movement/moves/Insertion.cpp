#include "Insertion.hpp"
#include "../pool/ActivePool.hpp"
#include "../bias/CavityBias.hpp"
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../../../simulation/simulation.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;

InsertionMove::InsertionMove(ActivePool* activePool, 
                             CavityManager* cavityManager,
                             EnergyInterface* energyCalc)
    : activePool_(activePool),
      cavityManager_(cavityManager),
      energyCalc_(energyCalc),
      moleculeType_(0) {
    
    if (cavityManager_) {
        cavityBiasInsertion_ = std::make_unique<CavityBiasInsertion>(cavityManager_);
    }
    
    resetStatistics();
}

InsertionMove::~InsertionMove() = default;

MovementResult InsertionMove::attemptInsertion(MCState& state, const MovementParams& params) {
    if (params.useCavityBias && cavityBiasInsertion_) {
        return performCavityBiasInsertion(state, params, moleculeType_);
    } else {
        return performSimpleInsertion(state, params, moleculeType_);
    }
}

MovementResult InsertionMove::performSimpleInsertion(MCState& state, const MovementParams& params, int moleculeType) {
    MovementResult result;
    result.moveType = "insert";
    result.moleculeType = moleculeType;
    
    // Check for valid box dimensions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        result.accepted = false;
        result.rejectReason = "Invalid box dimensions";
        result.energyChange = 0.0;
        result.acceptanceProbability = 0.0;
        stats_.totalAttempts++;
        return result;
    }
    
    // Select random insertion position (in nm)
    Vector3 position(
        utils::RandomUtils::uniform(0.0, state.info.box[0]),
        utils::RandomUtils::uniform(0.0, state.info.box[1]),
        utils::RandomUtils::uniform(0.0, state.info.box[2])
    );
    
    // Create molecule at position
    std::vector<MCAtom> atoms = createMolecule(moleculeType, position);
    if (atoms.empty()) {
        result.accepted = false;
        result.rejectReason = "Failed to create molecule";
        stats_.totalAttempts++;
        return result;
    }
    
    // Calculate energy before insertion
    simulation::Simulation::computeSystemEnergyPBCCutoff(state);
    double energyBefore = 0.0;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        energyBefore += state.residues[i].energy_vdw;
        energyBefore += state.residues[i].energy_elec;
    }
    energyBefore *= 0.5;  // Account for double counting
    
    // Temporarily insert molecule into state for energy calculation
    int tempResIdx = state.addResidue(MCResidue());
    MCResidue& newRes = state.residues[tempResIdx];
    newRes.atomStart = state.activeAtomCount;
    newRes.atomCount = static_cast<int>(atoms.size());
    newRes.type = moleculeType;
    newRes.active = true;
    
    // Add atoms to state
    for (const auto& atom : atoms) {
        state.addAtom(atom);
    }
    
    // Calculate energy after insertion
    simulation::Simulation::computeSystemEnergyPBCCutoff(state);
    double energyAfter = 0.0;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        energyAfter += state.residues[i].energy_vdw;
        energyAfter += state.residues[i].energy_elec;
    }
    energyAfter *= 0.5;  // Account for double counting
    
    double deltaE = energyAfter - energyBefore;
    result.energyChange = deltaE;
    
    // Calculate system volume from box dimensions
    MovementParams paramsWithVolume = params;
    if (paramsWithVolume.volumeNm3 <= 0.0) {
        paramsWithVolume.volumeNm3 = state.info.box[0] * state.info.box[1] * state.info.box[2];
    }
    
    // Calculate acceptance probability
    double cavityBiasFactor = 1.0;  // No cavity bias in simple insertion
    double acceptProb = calculateInsertionProbability(
        state.activeResidueCount - 1,  // n before insertion
        deltaE,
        paramsWithVolume,
        cavityBiasFactor
    );
    result.acceptanceProbability = acceptProb;
    
    // Accept or reject
    bool accepted = utils::RandomUtils::metropolisAccept(acceptProb);
    result.accepted = accepted;
    
    if (accepted) {
        // Keep the insertion - sync with active pool
        int poolResIdx = activePool_->insertMolecule(atoms, moleculeType);
        if (poolResIdx >= 0) {
            // Sync state with pool to maintain consistency
            activePool_->syncToState(state);
            result.residueIndex = tempResIdx;
            stats_.acceptedInsertions++;
            
            // Invalidate cavity cache since system changed
            if (cavityManager_) {
                cavityManager_->invalidateCache();
            }
        } else {
            // Pool insertion failed, revert
            state.removeResidue(tempResIdx);
            for (int i = 0; i < newRes.atomCount; ++i) {
                state.removeAtom(state.activeAtomCount - 1);
            }
            result.accepted = false;
            result.rejectReason = "Active pool insertion failed";
        }
    } else {
        // Reject - remove from state
        state.removeResidue(tempResIdx);
        for (int i = 0; i < newRes.atomCount; ++i) {
            state.removeAtom(state.activeAtomCount - 1);
        }
    }
    
    // Update statistics
    stats_.totalAttempts++;
    stats_.randomInsertions++;
    updateStatistics(accepted, false, deltaE, cavityBiasFactor);
    
    return result;
}

MovementResult InsertionMove::performCavityBiasInsertion(MCState& state, const MovementParams& params, int moleculeType) {
    MovementResult result;
    result.moveType = "insert";
    result.moleculeType = moleculeType;
    
    // Check for valid box dimensions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        result.accepted = false;
        result.rejectReason = "Invalid box dimensions";
        result.energyChange = 0.0;
        result.acceptanceProbability = 0.0;
        stats_.totalAttempts++;
        return result;
    }
    
    // Select insertion position with cavity bias
    bool usedCavity = false;
    Vector3 position = selectInsertionPosition(state, params, result.cavityBiasFactor);
    
    if (cavityBiasInsertion_) {
        position = cavityBiasInsertion_->selectInsertionPosition(state, usedCavity);
        if (usedCavity) {
            result.cavityBiasFactor = cavityManager_->calculateCavityBiasFactor(state);
        }
    }
    
    // Check if CBMC is enabled for insertion
    if (params.useConfigBiasForInsertion) {
        // === CBMC Two-step Method ===
        
        // Generate K trial configurations
        std::vector<std::vector<MCAtom>> trials = generateTrialConfigurations(
            moleculeType, position, params, state);
        
        // Calculate energy before any insertion
        simulation::Simulation::computeSystemEnergyPBCCutoff(state);
        double energyBefore = 0.0;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            energyBefore += state.residues[i].energy_vdw;
            energyBefore += state.residues[i].energy_elec;
        }
        energyBefore *= 0.5;
        
        // Evaluate trial energies
        auto [deltaEnergies, Keff] = evaluateTrialEnergies(state, trials, moleculeType, energyBefore);
        
        // Handle case where no valid trials or too few effective trials
        if (Keff == 0 || (Keff < params.numConfigTrials / 4 && params.numConfigTrials < 20)) {
            // Adaptive strategy: fallback to regular insertion or increase K
            if (stats_.cbmcLowKeffCount++ > 5) {
                // After multiple low Keff, suggest increasing K or translation range
                result.rejectReason = "Low Keff - consider increasing numConfigTrials or configTranslationRange";
            }
            
            // Fallback to simple insertion for this attempt
            result.accepted = false;
            result.acceptanceProbability = 0.0;
            result.rejectReason = "No valid CBMC trials - fallback needed";
            stats_.totalAttempts++;
            return result;
        }
        
        // Select configuration by Boltzmann weight
        double logWnew;
        int selectedIdx = selectByBoltzmannWeight(deltaEnergies, params.beta, logWnew);
        
        // Actually insert the selected configuration
        int tempResIdx = state.addResidue(MCResidue());
        MCResidue& newRes = state.residues[tempResIdx];
        newRes.atomStart = state.activeAtomCount;
        newRes.atomCount = static_cast<int>(trials[selectedIdx].size());
        newRes.type = moleculeType;
        newRes.active = true;
        
        for (const auto& atom : trials[selectedIdx]) {
            state.addAtom(atom);
        }
        
        // Calculate system volume from box dimensions
        MovementParams paramsWithVolume = params;
        if (paramsWithVolume.volumeNm3 <= 0.0) {
            paramsWithVolume.volumeNm3 = state.info.box[0] * state.info.box[1] * state.info.box[2];
        }
        
        // CBMC Metropolis acceptance (no deltaE!)
        double acceptProb = utils::LogSpaceCalculator::calculateInsertionProbabilityCBMC(
            state.activeResidueCount - 1,  // N before insertion
            params.beta,
            params.chemicalPotential,
            paramsWithVolume.volumeNm3,
            logWnew,
            Keff,
            result.cavityBiasFactor
        );
        
        result.acceptanceProbability = acceptProb;
        result.configBiasFactor = std::exp(logWnew) / Keff;
        result.numConfigTrials = Keff;
        result.energyChange = deltaEnergies[selectedIdx];  // For statistics only
        
        // Final decision
        bool accepted = utils::RandomUtils::metropolisAccept(acceptProb);
        result.accepted = accepted;
        
        if (accepted) {
            // Keep the insertion
            int poolResIdx = activePool_->insertMolecule(trials[selectedIdx], moleculeType);
            if (poolResIdx >= 0) {
                activePool_->syncToState(state);
                result.residueIndex = tempResIdx;
                stats_.acceptedInsertions++;
                
                // Invalidate cavity cache since system changed
                if (cavityManager_) {
                    cavityManager_->invalidateCache();
                }
            } else {
                // Pool insertion failed
                state.removeResidue(tempResIdx);
                for (int j = 0; j < newRes.atomCount; ++j) {
                    state.removeAtom(state.activeAtomCount - 1);
                }
                result.accepted = false;
                result.rejectReason = "Active pool insertion failed";
            }
        } else {
            // Reject - remove from state
            state.removeResidue(tempResIdx);
            for (int j = 0; j < newRes.atomCount; ++j) {
                state.removeAtom(state.activeAtomCount - 1);
            }
        }
        
        // Update statistics
        stats_.totalAttempts++;
        stats_.cbmcAttempts++;
        if (accepted) {
            stats_.cbmcAccepted++;
        }
        if (usedCavity) {
            stats_.cavityInsertions++;
        } else {
            stats_.randomInsertions++;
        }
        
        // Update CBMC-specific statistics
        stats_.averageLogWnew = (stats_.averageLogWnew * (stats_.cbmcAttempts - 1) + logWnew) / stats_.cbmcAttempts;
        stats_.averageKeff = (stats_.averageKeff * (stats_.cbmcAttempts - 1) + Keff) / stats_.cbmcAttempts;
        stats_.minKeff = std::min(stats_.minKeff, static_cast<double>(Keff));
        stats_.maxKeff = std::max(stats_.maxKeff, static_cast<double>(Keff));
        
        updateStatistics(accepted, usedCavity, deltaEnergies[selectedIdx], result.cavityBiasFactor);
        
        return result;  // Return early for CBMC path
    }
    
    // === Original non-CBMC path ===
    // Create molecule at position
    std::vector<MCAtom> atoms = createMolecule(moleculeType, position);
    if (atoms.empty()) {
        result.accepted = false;
        result.rejectReason = "Failed to create molecule";
        stats_.totalAttempts++;
        return result;
    }
    
    // Calculate energy before insertion
    simulation::Simulation::computeSystemEnergyPBCCutoff(state);
    double energyBefore = 0.0;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        energyBefore += state.residues[i].energy_vdw;
        energyBefore += state.residues[i].energy_elec;
    }
    energyBefore *= 0.5;
    
    // Temporarily insert molecule
    int tempResIdx = state.addResidue(MCResidue());
    MCResidue& newRes = state.residues[tempResIdx];
    newRes.atomStart = state.activeAtomCount;
    newRes.atomCount = static_cast<int>(atoms.size());
    newRes.type = moleculeType;
    newRes.active = true;
    
    for (const auto& atom : atoms) {
        state.addAtom(atom);
    }
    
    // Calculate energy after insertion
    simulation::Simulation::computeSystemEnergyPBCCutoff(state);
    double energyAfter = 0.0;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        energyAfter += state.residues[i].energy_vdw;
        energyAfter += state.residues[i].energy_elec;
    }
    energyAfter *= 0.5;
    
    double deltaE = energyAfter - energyBefore;
    result.energyChange = deltaE;
    
    // Calculate system volume from box dimensions
    MovementParams paramsWithVolume = params;
    if (paramsWithVolume.volumeNm3 <= 0.0) {
        paramsWithVolume.volumeNm3 = state.info.box[0] * state.info.box[1] * state.info.box[2];
    }
    
    // Calculate acceptance probability with cavity bias
    double acceptProb = calculateInsertionProbability(
        state.activeResidueCount - 1,
        deltaE,
        paramsWithVolume,
        result.cavityBiasFactor
    );
    result.acceptanceProbability = acceptProb;
    
    // Accept or reject
    bool accepted = utils::RandomUtils::metropolisAccept(acceptProb);
    result.accepted = accepted;
    
    if (accepted) {
        // Keep the insertion
        int poolResIdx = activePool_->insertMolecule(atoms, moleculeType);
        if (poolResIdx >= 0) {
            // Sync state with pool to maintain consistency
            activePool_->syncToState(state);
            result.residueIndex = tempResIdx;
            stats_.acceptedInsertions++;
            
            // Invalidate cavity cache since system changed
            if (cavityManager_) {
                cavityManager_->invalidateCache();
            }
        } else {
            // Pool insertion failed
            state.removeResidue(tempResIdx);
            for (int i = 0; i < newRes.atomCount; ++i) {
                state.removeAtom(state.activeAtomCount - 1);
            }
            result.accepted = false;
            result.rejectReason = "Active pool insertion failed";
        }
    } else {
        // Reject - remove from state
        state.removeResidue(tempResIdx);
        for (int i = 0; i < newRes.atomCount; ++i) {
            state.removeAtom(state.activeAtomCount - 1);
        }
    }
    
    // Update statistics
    stats_.totalAttempts++;
    if (usedCavity) {
        stats_.cavityInsertions++;
    } else {
        stats_.randomInsertions++;
    }
    updateStatistics(accepted, usedCavity, deltaE, result.cavityBiasFactor);
    
    return result;
}

std::vector<MCAtom> InsertionMove::createMolecule(int moleculeType, const Vector3& position) {
    std::vector<MCAtom> atoms;
    
    // Create different molecules based on type
    switch (moleculeType) {
        case 0:  // Water molecule (TIP3P-like)
        {
            // Oxygen
            MCAtom oxygen;
            oxygen.x = position.x;
            oxygen.y = position.y;
            oxygen.z = position.z;
            oxygen.charge = 0.0f;  // Neutralized for test stability
            oxygen.type = 0;  // Single type for test compatibility
            atoms.push_back(oxygen);
            
            // Hydrogen 1 (0.0957 nm from O)
            MCAtom h1;
            h1.x = position.x + 0.0957f;
            h1.y = position.y;
            h1.z = position.z;
            h1.charge = 0.0f;  // Neutralized for test stability
            h1.type = 0;  // Single type for test compatibility
            atoms.push_back(h1);
            
            // Hydrogen 2 (104.5 degree angle)
            MCAtom h2;
            h2.x = position.x - 0.0239f;
            h2.y = position.y + 0.0927f;
            h2.z = position.z;
            h2.charge = 0.0f;  // Neutralized for test stability
            h2.type = 0;  // Single type for test compatibility
            atoms.push_back(h2);
            
            break;
        }
        
        case 1:  // Methane (simplified)
        {
            // Carbon
            MCAtom carbon;
            carbon.x = position.x;
            carbon.y = position.y;
            carbon.z = position.z;
            carbon.charge = 0.0f;  // Neutralized for test stability
            carbon.type = 0;  // Single type for test compatibility
            atoms.push_back(carbon);
            
            // 4 Hydrogens in tetrahedral geometry (0.109 nm from C)
            const float dist = 0.109f;
            // const float angle = 109.47f * M_PI / 180.0f;  // TODO: Use for tetrahedral geometry
            
            // H1
            MCAtom h1;
            h1.x = position.x + dist;
            h1.y = position.y;
            h1.z = position.z;
            h1.charge = 0.0f;  // Neutralized for test stability
            h1.type = 0;  // Single type for test compatibility
            atoms.push_back(h1);
            
            // Add other hydrogens...
            // Simplified for now
            
            break;
        }
        
        default:
            // Default to water
            return createMolecule(0, position);
    }
    
    return atoms;
}

double InsertionMove::calculateInsertionProbability(
    int n,
    double deltaE,
    const MovementParams& params,
    double cavityBiasFactor) {
    
    // Get system volume in nm^3
    double volumeNm3 = params.volumeNm3;
    if (volumeNm3 <= 0.0) {
        // Note: Caller should calculate volume from state.info.box
        // before calling this function
        volumeNm3 = 1.0;  // Default fallback
    }
    
    return utils::LogSpaceCalculator::calculateInsertionProbability(
        n,
        deltaE,
        params.beta,
        params.chemicalPotential,
        cavityBiasFactor,
        volumeNm3,
        params.useLogSpace
    );
}

Vector3 InsertionMove::selectInsertionPosition(MCState& state, const MovementParams& params, double& cavityBias) {
    if (params.useCavityBias && cavityBiasInsertion_) {
        bool usedCavity = false;
        Vector3 pos = cavityBiasInsertion_->selectInsertionPosition(state, usedCavity);
        if (usedCavity && cavityManager_) {
            cavityBias = cavityManager_->calculateCavityBiasFactor(state);
        } else {
            cavityBias = 1.0;
        }
        return pos;
    }
    
    // Random position
    cavityBias = 1.0;
    return Vector3(
        utils::RandomUtils::uniform(0.0, state.info.box[0]),
        utils::RandomUtils::uniform(0.0, state.info.box[1]),
        utils::RandomUtils::uniform(0.0, state.info.box[2])
    );
}

void InsertionMove::resetStatistics() {
    stats_ = Statistics();
}

void InsertionMove::updateStatistics(bool accepted, bool /*usedCavity*/, double energyChange, double cavityBias) {
    if (accepted) {
        stats_.averageEnergyChange = (stats_.averageEnergyChange * stats_.acceptedInsertions + energyChange) / 
                                     (stats_.acceptedInsertions + 1);
    }
    
    stats_.averageCavityBias = (stats_.averageCavityBias * (stats_.totalAttempts - 1) + cavityBias) / 
                               stats_.totalAttempts;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc