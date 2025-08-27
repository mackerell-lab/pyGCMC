#include "MovementMain.hpp"
#include "../pool/ActivePool.hpp"
#include "../bias/CavityBias.hpp"
#include "../bias/ConfigBias.hpp"
#include "../moves/Insertion.hpp"
#include "../moves/Deletion.hpp"
#include "../moves/Translation.hpp"
#include "../moves/Rotation.hpp"
#include "../moves/MultiInsertionCBMC.hpp"
#include "../common/MovementUtils.hpp"
#ifdef PYGCMC_USE_PROPOSAL_LAYER
#include "../proposal/ProposalMain.hpp"
#endif
#include <random>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Private implementation class
class MovementModule::Impl {
public:
    std::unique_ptr<InsertionMove> insertionMove;
    std::unique_ptr<DeletionMove> deletionMove;
    std::unique_ptr<TranslationMove> translationMove;
    std::unique_ptr<RotationMove> rotationMove;
    std::unique_ptr<MultiInsertionCBMC> multiInsertionCBMC;
    
    std::mt19937 rng;
    
    Impl() : rng(std::chrono::steady_clock::now().time_since_epoch().count()) {}
};

MovementModule::MovementModule()
    : MovementModule(MovementParams()) {
    // Delegate to parameterized constructor with default params
}

MovementModule::MovementModule(const MovementParams& params)
    : pImpl_(std::make_unique<Impl>()),
      params_(params),
      activePool_(std::make_unique<ActivePool>(params.maxAtoms, params.maxResidues)),
      cavityManager_(std::make_unique<CavityManager>(params.cavityGridSpacing * 10.0, params.probeRadius * 10.0)),
      configBiasManager_(std::make_unique<ConfigBiasManager>(params.numConfigTrials)),
      energyCalc_(nullptr) {
    
    initializeComponents();
}

MovementModule::~MovementModule() = default;

void MovementModule::initializeComponents() {
    // Set random seed if specified
    if (params_.seed != 0) {
        utils::RandomUtils::setSeed(params_.seed);
    }
    
    // Adjust cavity grid spacing to align with independence scale
    if (cavityManager_ && params_.useCavityBias) {
        double spacingNm = std::max(params_.cavityGridSpacing, params_.minRegionSeparationNm * 0.5);
        // Avoid overly fine grids that cause performance issues
        spacingNm = std::max(spacingNm, 0.2);  // Minimum 0.2 nm
        cavityManager_->setGridSpacing(spacingNm * 10.0);  // Convert to Angstroms
    }
    
    // Initialize movement implementations
    pImpl_->insertionMove = std::make_unique<InsertionMove>(
        activePool_.get(), cavityManager_.get(), energyCalc_.get());
    
    pImpl_->deletionMove = std::make_unique<DeletionMove>(
        activePool_.get(), energyCalc_.get());
    
    pImpl_->translationMove = std::make_unique<TranslationMove>(
        activePool_.get(), energyCalc_.get());
    
    pImpl_->rotationMove = std::make_unique<RotationMove>(
        activePool_.get(), configBiasManager_.get(), energyCalc_.get());
    
    // Initialize multi-insertion CBMC if enabled
    if (params_.useMultiInsertionCBMC && params_.maxParallelInsertions > 0) {
        MultiInsertionConfig config;
        config.numTrialsPerRegion = params_.numConfigTrials;
        config.maxParallelInsertions = params_.maxParallelInsertions;
        config.minSeparation = params_.minRegionSeparationNm;  // Keep in nm
        config.chemicalPotential = params_.chemicalPotential;
        config.displacementFraction = params_.multiDisplacementFraction;
        config.useRegionVolume = params_.multiUseRegionVolume;
        config.useCavityBias = params_.useCavityBias;  // Pass cavity bias flag
        
        // Set independence control parameters
        // Note: cutoff will be set from MCState when actually used
        config.cutoffNm = 1.2;  // Default cutoff in nm (will be updated from state)
        config.moleculeExtentNm = 0.15;    // Default for water, should be configurable
        config.enforceIndependence = true; // Always enforce for correctness
        config.recomputeAfterAccept = false; // Expensive fallback, off by default
        
        pImpl_->multiInsertionCBMC = std::make_unique<MultiInsertionCBMC>(config, cavityManager_.get());
        if (params_.seed != 0) {
            pImpl_->multiInsertionCBMC->setSeed(params_.seed);
        }
    }
    
    // Initialize statistics
    stats_["insert"] = Statistics();
    stats_["delete"] = Statistics();
    stats_["translate"] = Statistics();
    stats_["rotate"] = Statistics();
    stats_["multi_insert"] = Statistics();
}

MovementResult MovementModule::attemptInsertion(MCState& state, int moleculeType) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    MovementResult result;
    
    if (params_.useCavityBias) {
        result = pImpl_->insertionMove->performCavityBiasInsertion(state, params_, moleculeType);
    } else {
        result = pImpl_->insertionMove->performSimpleInsertion(state, params_, moleculeType);
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    result.computeTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    updateStatistics("insert", result.accepted, result.energyChange);
    
    // Remember last accepted insertion for paired deletion
    if (result.accepted) {
        lastInsertedResidueIndex_ = result.residueIndex;
    }
    
    return result;
}

MovementResult MovementModule::attemptDeletion(MCState& state, int residueIndex) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Prefer deleting the last accepted insertion when user didn't specify an index
    if (residueIndex < 0 && lastInsertedResidueIndex_ >= 0 &&
        lastInsertedResidueIndex_ < state.activeResidueCount &&
        state.residues[lastInsertedResidueIndex_].active) {
        residueIndex = lastInsertedResidueIndex_;
        lastInsertedResidueIndex_ = -1; // consume once
    }
    
    MovementResult result = pImpl_->deletionMove->performDeletion(state, params_, residueIndex);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    result.computeTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    updateStatistics("delete", result.accepted, result.energyChange);
    
    // Invalidate cavity cache after accepted deletion
    if (result.accepted && cavityManager_) {
        cavityManager_->invalidateCache();
    }
    
    return result;
}

MovementResult MovementModule::attemptTranslation(MCState& state, int residueIndex) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    MovementResult result = pImpl_->translationMove->performTranslation(state, params_, residueIndex);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    result.computeTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    updateStatistics("translate", result.accepted, result.energyChange);
    
    // Invalidate cavity cache after accepted translation
    if (result.accepted && cavityManager_) {
        cavityManager_->invalidateCache();
    }
    
    return result;
}

MovementResult MovementModule::attemptRotation(MCState& state, int residueIndex) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    MovementResult result;
    
    if (params_.useConfigBias) {
        result = pImpl_->rotationMove->performConfigBiasRotation(state, params_, residueIndex);
    } else {
        result = pImpl_->rotationMove->performSimpleRotation(state, params_, residueIndex);
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    result.computeTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    updateStatistics("rotate", result.accepted, result.energyChange);
    
    // Invalidate cavity cache after accepted rotation
    if (result.accepted && cavityManager_) {
        cavityManager_->invalidateCache();
    }
    
    return result;
}

MovementResult MovementModule::attemptCavityBiasInsertion(MCState& state, int moleculeType) {
    // Force cavity bias
    bool originalSetting = params_.useCavityBias;
    params_.useCavityBias = true;
    
    MovementResult result = attemptInsertion(state, moleculeType);
    
    params_.useCavityBias = originalSetting;
    
    return result;
}

MovementResult MovementModule::attemptConfigBiasRotation(MCState& state, int residueIndex) {
    // Force configurational bias
    bool originalSetting = params_.useConfigBias;
    params_.useConfigBias = true;
    
    MovementResult result = attemptRotation(state, residueIndex);
    
    params_.useConfigBias = originalSetting;
    
    return result;
}

std::vector<Vector3> MovementModule::findCavities(const MCState& state) {
    return cavityManager_->findCavities(state);
}

double MovementModule::calculateAcceptanceRate(const std::string& moveType) const {
    auto it = stats_.find(moveType);
    if (it != stats_.end()) {
        return it->second.acceptanceRate();
    }
    return 0.0;
}

void MovementModule::resetStatistics() {
    for (auto& pair : stats_) {
        pair.second = Statistics();
    }
    
    pImpl_->insertionMove->resetStatistics();
    pImpl_->deletionMove->resetStatistics();
    pImpl_->translationMove->resetStatistics();
    pImpl_->rotationMove->resetStatistics();
    
    activePool_->resetStatistics();
    cavityManager_->resetStatistics();
    configBiasManager_->resetStatistics();
}

void MovementModule::setParams(const MovementParams& params) {
    params_ = params;
    // Ensure derived values and validation are up-to-date
    params_.updateDerivedParameters();
    
    // Reseed RNGs when a non-zero seed is provided
    if (params_.seed != 0) {
        utils::RandomUtils::setSeed(params_.seed);
        if (pImpl_ && pImpl_->multiInsertionCBMC) {
            pImpl_->multiInsertionCBMC->setSeed(params_.seed);
        }
    }
    
    // Update component configurations (convert nm to Angstroms by multiplying by 10)
    cavityManager_->setGridSpacing(params_.cavityGridSpacing * 10.0);
    cavityManager_->setProbeRadius(params_.probeRadius * 10.0);
    configBiasManager_->setNumTrials(params_.numConfigTrials);
    configBiasManager_->setTranslationRange(params_.configTranslationRange * 10.0);
    activePool_->setFragmentationThreshold(params_.fragmentationThreshold);
}

MovementParams MovementModule::getParams() const {
    return params_;
}

std::map<std::string, MovementModule::Statistics> MovementModule::getStatistics() const {
    return stats_;
}

std::vector<MovementResult> MovementModule::attemptMultiInsertionCBMC(MCState& state, int moleculeType) {
    std::vector<MovementResult> results;
    
    if (!pImpl_->multiInsertionCBMC) {
        // Initialize if not already done
        MultiInsertionConfig config;
        config.numTrialsPerRegion = params_.numConfigTrials > 0 ? params_.numConfigTrials : 10;
        // Ensure at least 1 parallel insertion for fallback path
        config.maxParallelInsertions = std::max(1, params_.maxParallelInsertions);
        config.minSeparation = params_.minRegionSeparationNm > 0 ? params_.minRegionSeparationNm : 1.5;  // Default 1.5 nm
        config.chemicalPotential = params_.chemicalPotential;
        config.displacementFraction = params_.multiDisplacementFraction > 0 ? params_.multiDisplacementFraction : 0.5;
        config.useRegionVolume = params_.multiUseRegionVolume;
        config.useCavityBias = params_.useCavityBias;  // Pass cavity bias flag
        
        // Set independence control parameters
        // Get cutoff from state
        config.cutoffNm = state.info.cutoff;  // Already in nm
        config.moleculeExtentNm = 0.15;    // Default for water
        config.enforceIndependence = true; // Always enforce for correctness
        config.recomputeAfterAccept = false; // Expensive fallback, off by default
        
        pImpl_->multiInsertionCBMC = std::make_unique<MultiInsertionCBMC>(config, cavityManager_.get());
        
        // Set seed for reproducibility in fallback path
        if (params_.seed != 0) {
            pImpl_->multiInsertionCBMC->setSeed(params_.seed);
        }
    }
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Perform multi-insertion
    auto [acceptCount, regions] = pImpl_->multiInsertionCBMC->performMultiInsertion(
        state, moleculeType, params_);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    double totalTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    // Create results for each region
    // First, commit accepted molecules to the state in batch
    for (const auto& region : regions) {
        if (region.accepted && region.selectedConfig >= 0) {
            // Insert accepted configuration into active pool and state
            activePool_->insertMolecule(region.trialConfigs[region.selectedConfig], moleculeType);
            
            // Sync to state - add atoms and residue
            for (const auto& atom : region.trialConfigs[region.selectedConfig]) {
                state.addAtom(atom);
            }
            
            // Add residue
            MCResidue newRes;
            newRes.atomStart = state.activeAtomCount - region.trialConfigs[region.selectedConfig].size();
            newRes.atomCount = region.trialConfigs[region.selectedConfig].size();
            newRes.type = moleculeType;
            newRes.active = true;
            state.addResidue(newRes);
        }
    }
    
    // Count accepts
    int actualAccepts = 0;
    for (const auto& region : regions) {
        if (region.accepted) actualAccepts++;
    }
    
    // Invalidate cavity cache only when necessary (when molecules were actually accepted)
    if (cavityManager_ && params_.useCavityBias && actualAccepts > 0) {
        cavityManager_->invalidateCache();
    }
    
    // Create result entries
    for (const auto& region : regions) {
        MovementResult result;
        result.moveType = "multi_insert";
        result.accepted = region.accepted;
        result.energyChange = region.accepted ? region.trialEnergies[region.selectedConfig] : 0.0;
        result.residueIndex = -1;  // Multiple residues
        result.computeTimeMs = totalTimeMs / regions.size();  // Average time per region
        results.push_back(result);
        
        updateStatistics("multi_insert", result.accepted, result.energyChange);
    }
    
    return results;
}

void MovementModule::updateStatistics(const std::string& moveType, bool accepted, double energyChange) {
    auto& stat = stats_[moveType];
    stat.attempts++;
    if (accepted) {
        stat.accepts++;
    }
    stat.totalEnergyChange += energyChange;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc