#include "GCMCModule.hpp"
#include "GCMCEngine.hpp"
#include "GCMCMoveSelector.hpp"
#include "GCMCBias.hpp"
#include "GCMCStats.hpp"
#include "GCMCAcceptance.hpp"
#include "../reservoir/fragment_reservoir.hpp"
#include "../bias/CavityBias.hpp"
#include "../bias/ConfigBias.hpp"
#include "../../energy/EnergyModule.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Normalize move probabilities to sum to 1.0
void GCMCModule::Config::normalizeProbs() {
    double sum = insertProb + deleteProb + translateProb + rotateProb + swapProb;
    if (sum > 0) {
        insertProb /= sum;
        deleteProb /= sum;
        translateProb /= sum;
        rotateProb /= sum;
        swapProb /= sum;
    }
}

// Validate configuration
void GCMCModule::Config::validate() const {
    // Temperature validation
    if (temperature <= 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): temperature must be > 0, got " 
                                   + std::to_string(temperature));
    }
    
    // Pressure validation
    if (pressure <= 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): pressure must be > 0, got " 
                                   + std::to_string(pressure));
    }
    
    // Step validation
    if (equilibrationSteps < 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): equilibrationSteps must be >= 0, got " 
                                   + std::to_string(equilibrationSteps));
    }
    if (productionSteps <= 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): productionSteps must be > 0, got " 
                                   + std::to_string(productionSteps));
    }
    if (saveFrequency < 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): saveFrequency must be >= 0, got " 
                                   + std::to_string(saveFrequency));
    }
    
    // Probability validation
    double probSum = insertProb + deleteProb + translateProb + rotateProb + swapProb;
    if (probSum <= 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): sum of move probabilities must be > 0, got " 
                                   + std::to_string(probSum));
    }
    if (insertProb < 0 || deleteProb < 0 || translateProb < 0 || rotateProb < 0 || swapProb < 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): all move probabilities must be >= 0");
    }
    
    // Bias parameter validation
    if (gridSpacing <= 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): gridSpacing must be > 0, got " 
                                   + std::to_string(gridSpacing));
    }
    if (probeRadius < 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): probeRadius must be >= 0, got " 
                                   + std::to_string(probeRadius));
    }
    if (configTrials < 1) {
        throw std::invalid_argument("GCMCModule::Config::validate(): configTrials must be >= 1, got " 
                                   + std::to_string(configTrials));
    }
    
    // Cutoff validation
    if (cutoff <= 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): cutoff must be > 0, got " 
                                   + std::to_string(cutoff));
    }
    
    // Cluster cutoff validation
    if (useClusterMoves && clusterCutoff <= 0) {
        throw std::invalid_argument("GCMCModule::Config::validate(): clusterCutoff must be > 0 when useClusterMoves is true, got " 
                                   + std::to_string(clusterCutoff));
    }
}

// Constructor
GCMCModule::GCMCModule(const Config& config) 
    : state_(nullptr),
      config_(config),
      initialized_(false),
      currentStep_(0),
      currentPhase_(Phase::EQUILIBRATION) {
    
    // Validate configuration
    config_.validate();
    
    // Normalize probabilities
    config_.normalizeProbs();
    
    // Create core components
    engine_ = std::make_unique<GCMCEngine>();
    
    // Configure reservoir to prevent slot reordering/compaction
    FragmentReservoir::Config reservoirConfig;
    reservoirConfig.autoCompact = false;          // CRITICAL: Never compact/reorder slots
    reservoirConfig.ghostRecycleRatio = 1.0;      // Always recycle ghosts in-place
    reservoirConfig.maxGhosts = 10000;            // Large limit to avoid purging
    reservoirConfig.maxInstances = 10000;         // Sufficient for most simulations
    reservoir_ = std::make_unique<FragmentReservoir>(reservoirConfig);
    
    moveSelector_ = std::make_unique<GCMCMoveSelector>();
    biasCalc_ = std::make_unique<GCMCBias>();
    acceptCalc_ = std::make_unique<GCMCAcceptance>();
    statistics_ = std::make_unique<GCMCStats>();
    
    // Initialize cavity manager if needed
    if (config_.useCavityBias) {
        cavityManager_ = std::make_unique<CavityManager>();
    }
    
    // Initialize config bias if needed
    if (config_.useConfigBias) {
        configBias_ = std::make_unique<ConfigBiasManager>();
    }
    
    // Configure move selector
    std::map<GCMCMoveSelector::MoveType, double> probs;
    probs[GCMCMoveSelector::MoveType::INSERT] = config_.insertProb;
    probs[GCMCMoveSelector::MoveType::DELETE] = config_.deleteProb;
    probs[GCMCMoveSelector::MoveType::TRANSLATE] = config_.translateProb;
    probs[GCMCMoveSelector::MoveType::ROTATE] = config_.rotateProb;
    probs[GCMCMoveSelector::MoveType::SWAP] = config_.swapProb;
    moveSelector_->setProbabilities(probs);
    
    // Configure acceptance calculator
    acceptCalc_->setTemperature(config_.temperature);
}

// Destructor
GCMCModule::~GCMCModule() = default;

// Initialize with system state
void GCMCModule::initialize(model::montecarlo::MCState& state) {
    state_ = &state;
    
    // CRITICAL: Ensure periodicBox is initialized from info.box if not set
    // This prevents segfaults in functions that directly index periodicBox[0..2]
    if (state.periodicBox.size() < 3 && 
        state.info.box[0] > 0 && state.info.box[1] > 0 && state.info.box[2] > 0) {
        state.periodicBox.resize(3);
        state.periodicBox[0] = state.info.box[0];
        state.periodicBox[1] = state.info.box[1];
        state.periodicBox[2] = state.info.box[2];
        
        if (config_.verbose) {
            std::cout << "Auto-initialized periodicBox from info.box: " 
                     << state.periodicBox[0] << " x " 
                     << state.periodicBox[1] << " x " 
                     << state.periodicBox[2] << " nm" << std::endl;
        }
    }
    
    // Initialize engine
    engine_->initialize(state_, reservoir_.get());
    engine_->setTemperature(config_.temperature);
    
    // Determine actual energy method to use based on system setup
    EnergyMethod actualMethod = config_.energyMethod;
    bool hasValidBox = (state.periodicBox.size() == 3 && 
                       state.periodicBox[0] > 0 && 
                       state.periodicBox[1] > 0 && 
                       state.periodicBox[2] > 0);
    
    // PME/Ewald require valid periodic box - fallback to DIRECT if not available
    if ((actualMethod == EnergyMethod::PME || actualMethod == EnergyMethod::EWALD) && !hasValidBox) {
        if (config_.verbose) {
            std::cerr << "WARNING: " << (actualMethod == EnergyMethod::PME ? "PME" : "Ewald") 
                     << " energy method requires valid periodic box. Falling back to DIRECT." << std::endl;
        }
        actualMethod = EnergyMethod::DIRECT;
    }
    
    engine_->setEnergyMethod(actualMethod);
    engine_->setCutoff(config_.cutoff);
    
    // Configure energy callback
    if (auto* callback = engine_->getEnergyCallback()) {
        callback->setEnergyMethod(actualMethod);
        callback->setParameters(true, hasValidBox);  // Use cutoff, PBC only if box valid
        
        // Initialize Ewald/PME if needed and possible
        if (actualMethod == EnergyMethod::EWALD && hasValidBox) {
            callback->initializeEwald(config_.cutoff, state.periodicBox);
        } else if (actualMethod == EnergyMethod::PME && hasValidBox) {
            std::vector<int> gridSize = {64, 64, 64};  // Default PME grid
            callback->initializePME(config_.cutoff, state.periodicBox, gridSize);
        }
    }
    
    // Setup acceptance calculator for proper GCMC
    if (state.periodicBox.size() == 3) {
        double volume = state.periodicBox[0] * state.periodicBox[1] * state.periodicBox[2];
        acceptCalc_->setVolume(volume);
    }
    acceptCalc_->setTemperature(config_.temperature);
    engine_->setAcceptanceCalculator(acceptCalc_.get());
    
    // Set bias components in engine
    if (cavityManager_) {
        engine_->setCavityManager(cavityManager_.get());
        
        // Initialize cavity manager with box dimensions
        if (state.periodicBox.size() == 3) {
            // Configure cavity manager - convert nm to Angstrom
            cavityManager_->setGridSpacing(config_.gridSpacing * 10.0);  // nm to Å
            cavityManager_->setProbeRadius(config_.probeRadius * 10.0);   // nm to Å
            // findCavities will automatically initialize the grid with the correct box size
            cavityManager_->findCavities(state);
        }
    }
    
    if (configBias_) {
        engine_->setConfigBiasManager(configBias_.get());
    }
    
    // Initialize bias calculator
    biasCalc_->initialize(state_);
    biasCalc_->setTemperature(config_.temperature);
    if (cavityManager_) {
        biasCalc_->setCavityManager(cavityManager_.get());
    }
    if (configBias_) {
        biasCalc_->setConfigBiasManager(configBias_.get());
    }
    
    // Initialize energy calculation based on method
    if (config_.energyMethod == EnergyMethod::PME && state.periodicBox.size() == 3) {
        double box[3] = {
            state.periodicBox[0],
            state.periodicBox[1],
            state.periodicBox[2]
        };
        initializePMEParameters(config_.cutoff, box);
    } else if (config_.energyMethod == EnergyMethod::EWALD && state.periodicBox.size() == 3) {
        double box[3] = {
            state.periodicBox[0],
            state.periodicBox[1],
            state.periodicBox[2]
        };
        initializeEwaldParameters(config_.cutoff, box);
    }
    
    // Initialize statistics
    statistics_->initialize(reservoir_->getTemplateCount());
    
    initialized_ = true;
}

// Add fragment type
int GCMCModule::addFragmentType(const FragmentTemplate& tmpl) {
    int typeId = reservoir_->addTemplate(tmpl);
    
    // Update acceptance calculator
    acceptCalc_->setChemicalPotential(typeId, tmpl.chemicalPotential);
    acceptCalc_->setActivity(typeId, tmpl.activity);
    
    // Update statistics
    statistics_->setFragmentInfo(typeId, tmpl.name, tmpl.chemicalPotential);
    
    return typeId;
}

// Set chemical potential
void GCMCModule::setChemicalPotential(int typeId, double mu) {
    FragmentTemplate* tmpl = const_cast<FragmentTemplate*>(reservoir_->getTemplate(typeId));
    if (tmpl) {
        tmpl->chemicalPotential = mu;
        tmpl->calculateActivity(config_.temperature);
        acceptCalc_->setChemicalPotential(typeId, mu);
        acceptCalc_->setActivity(typeId, tmpl->activity);
    }
}

// Set activity
void GCMCModule::setActivity(int typeId, double activity) {
    FragmentTemplate* tmpl = const_cast<FragmentTemplate*>(reservoir_->getTemplate(typeId));
    if (tmpl) {
        tmpl->activity = activity;
        tmpl->chemicalPotential = 8.314e-3 * config_.temperature * std::log(activity);
        acceptCalc_->setActivity(typeId, activity);
        acceptCalc_->setChemicalPotential(typeId, tmpl->chemicalPotential);
    }
}

// Run equilibration
void GCMCModule::runEquilibration(int steps) {
    if (!initialized_) {
        throw std::runtime_error("GCMCModule not initialized");
    }
    
    currentPhase_ = Phase::EQUILIBRATION;
    int nSteps = (steps < 0) ? config_.equilibrationSteps : steps;
    
    if (config_.verbose) {
        std::cout << "Starting equilibration for " << nSteps << " steps..." << std::endl;
    }
    
    runSteps(nSteps);
    
    if (config_.verbose) {
        std::cout << "Equilibration complete." << std::endl;
        printStatistics();
    }
}

// Run production
void GCMCModule::runProduction(int steps) {
    if (!initialized_) {
        throw std::runtime_error("GCMCModule not initialized");
    }
    
    currentPhase_ = Phase::PRODUCTION;
    int nSteps = (steps < 0) ? config_.productionSteps : steps;
    
    if (config_.verbose) {
        std::cout << "Starting production for " << nSteps << " steps..." << std::endl;
    }
    
    // Reset statistics for production phase
    statistics_->resetMoveStatistics();
    
    runSteps(nSteps);
    
    // Calculate final averages
    statistics_->updateAverages();
    statistics_->calculateFluctuations();
    
    if (state_->periodicBox.size() == 3) {
        double volume = state_->periodicBox[0] * state_->periodicBox[1] * state_->periodicBox[2];
        statistics_->calculateChemicalPotentials(volume, config_.temperature);
    }
    
    if (config_.verbose) {
        std::cout << "Production complete." << std::endl;
        printStatistics();
    }
}

// Run specified number of steps
void GCMCModule::runSteps(int nSteps) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < nSteps; ++i) {
        currentStep_++;
        
        // Perform move
        performMove();
        
        // Update statistics
        if (currentPhase_ == Phase::PRODUCTION) {
            updateStatistics();
        }
        
        // Save trajectory if needed
        if (config_.saveFrequency > 0 && currentStep_ % config_.saveFrequency == 0) {
            saveSnapshot();
        }
        
        // Verbose output
        if (config_.verbose && currentStep_ % 1000 == 0) {
            std::cout << "Step " << currentStep_ 
                     << " | Molecules: " << reservoir_->getActiveCount()
                     << " | Energy: " << calculateSystemEnergy()
                     << " | Accept rate: " << statistics_->getTranslateStats().acceptanceRate()
                     << std::endl;
        }
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    statistics_->recordStepTime(duration.count() / static_cast<double>(nSteps));
}

// Set random seed
void GCMCModule::setSeed(unsigned int seed) {
    if (engine_) {
        engine_->setSeed(seed);
    }
    if (moveSelector_) {
        moveSelector_->setSeed(seed + 1);  // Use different seed for move selector
    }
    // Also set seed for any other RNG components if needed
}

// Perform a move
bool GCMCModule::performMove() {
    // Select move type
    GCMCMoveSelector::MoveType moveType = moveSelector_->selectMove();
    
    bool accepted = false;
    
    switch (moveType) {
        case GCMCMoveSelector::MoveType::INSERT:
            accepted = attemptInsertion();
            break;
        case GCMCMoveSelector::MoveType::DELETE:
            accepted = attemptDeletion();
            break;
        case GCMCMoveSelector::MoveType::TRANSLATE:
            accepted = attemptTranslation();
            break;
        case GCMCMoveSelector::MoveType::ROTATE:
            accepted = attemptRotation();
            break;
        case GCMCMoveSelector::MoveType::SWAP:
            accepted = attemptSwap();
            break;
        case GCMCMoveSelector::MoveType::REGROWTH:
            if (config_.useRegrowth) {
                accepted = attemptRegrowth();
            }
            break;
        case GCMCMoveSelector::MoveType::CLUSTER:
            if (config_.useClusterMoves) {
                accepted = attemptClusterMove();
            }
            break;
        default:
            break;
    }
    
    // Record move statistics
    moveSelector_->recordAttempt(moveType);
    if (accepted) {
        moveSelector_->recordAcceptance(moveType);
    }
    
#ifdef DEBUG
    // Verify state consistency after move
    validateStateConsistency();
#endif
    
    return accepted;
}

// Attempt insertion
bool GCMCModule::attemptInsertion(int typeId) {
    // Select fragment type if not specified
    if (typeId < 0) {
        typeId = engine_->selectRandomFragment();
        if (typeId < 0) return false;
    }
    
    // Perform insertion
    auto result = engine_->attemptInsertion(typeId);
    
    // Update statistics
    statistics_->recordInsertion(typeId, result.accepted, result.deltaE);
    
    // Update cavity manager if needed
    if (result.accepted && cavityManager_) {
        cavityManager_->updateAfterInsertion(result.residueIndex, *state_);
    }
    
    return result.accepted;
}

// Attempt deletion
bool GCMCModule::attemptDeletion(int typeId) {
    // Select fragment type if not specified
    if (typeId < 0) {
        typeId = engine_->selectRandomFragment();
        if (typeId < 0) return false;
    }
    
    // Select instance to delete
    int residueIdx = engine_->selectRandomInstance(typeId);
    if (residueIdx < 0) return false;
    
    // Get position before deletion for cavity update
    auto instance = reservoir_->getInstanceByResidueIndex(residueIdx);
    movement::Vector3 position = instance ? instance->position : movement::Vector3();
    
    // Perform deletion
    auto result = engine_->attemptDeletion(typeId);
    
    // Update statistics
    statistics_->recordDeletion(typeId, result.accepted, result.deltaE);
    
    // Update cavity manager if needed
    if (result.accepted && cavityManager_) {
        cavityManager_->updateAfterDeletion(position, *state_);
    }
    
    return result.accepted;
}

// Attempt translation
bool GCMCModule::attemptTranslation() {
    // Select random instance
    int residueIdx = engine_->selectRandomInstance();
    if (residueIdx < 0) return false;
    
    // Perform translation
    auto result = engine_->attemptTranslation(residueIdx);
    
    // Update statistics
    statistics_->recordTranslation(result.fragmentType, result.accepted, result.deltaE);
    
    return result.accepted;
}

// Attempt rotation
bool GCMCModule::attemptRotation() {
    // Select random instance
    int residueIdx = engine_->selectRandomInstance();
    if (residueIdx < 0) return false;
    
    // Perform rotation
    auto result = engine_->attemptRotation(residueIdx);
    
    // Update statistics
    statistics_->recordRotation(result.fragmentType, result.accepted, result.deltaE);
    
    return result.accepted;
}

// Attempt swap
bool GCMCModule::attemptSwap() {
    // Select two different fragment types
    if (reservoir_->getTemplateCount() < 2) return false;
    
    int type1 = engine_->selectRandomFragment();
    int type2 = engine_->selectRandomFragment();
    while (type2 == type1) {
        type2 = engine_->selectRandomFragment();
    }
    
    // Perform swap
    auto result = engine_->attemptSwap(type1, type2);
    
    // Update statistics
    statistics_->recordSwap(type1, type2, result.accepted, result.deltaE);
    
    return result.accepted;
}

// Attempt regrowth
bool GCMCModule::attemptRegrowth() {
    // Select random instance
    int residueIdx = engine_->selectRandomInstance();
    if (residueIdx < 0) return false;
    
    // Perform regrowth
    auto result = engine_->attemptRegrowth(residueIdx);
    
    // Update statistics (record as rotation for now)
    statistics_->recordRotation(result.fragmentType, result.accepted, result.deltaE);
    
    return result.accepted;
}

// Attempt cluster move
bool GCMCModule::attemptClusterMove() {
    // Select seed instance
    int seedIdx = engine_->selectRandomInstance();
    if (seedIdx < 0) return false;
    
    // Perform cluster move
    auto result = engine_->attemptClusterMove(seedIdx, config_.clusterCutoff);
    
    // Update statistics (record as translation for now)
    statistics_->recordTranslation(result.fragmentType, result.accepted, result.deltaE);
    
    return result.accepted;
}

// Print statistics
void GCMCModule::printStatistics() const {
    statistics_->print();
}

// Save statistics to file
void GCMCModule::saveStatistics(const std::string& filename) const {
    statistics_->saveToFile(filename);
}

// Save snapshot
void GCMCModule::saveSnapshot() {
    if (config_.trajectoryFile.empty()) return;
    
    std::ofstream file(config_.trajectoryFile, std::ios::app);
    if (!file.is_open()) return;
    
    file << "MODEL " << currentStep_ << "\n";
    
    // Write all active atoms
    int atomIdx = 1;
    for (const auto& residue : state_->residues) {
        if (!residue.active) continue;
        
        for (const auto& atom : residue.atoms) {
            file << "ATOM  " << std::setw(5) << atomIdx++ 
                 << " " << std::setw(4) << atom.name
                 << " " << std::setw(3) << residue.resname
                 << "  " << std::setw(4) << residue.resid
                 << "    "
                 << std::fixed << std::setprecision(3)
                 << std::setw(8) << atom.position.x * 10.0  // nm to Angstrom
                 << std::setw(8) << atom.position.y * 10.0
                 << std::setw(8) << atom.position.z * 10.0
                 << "\n";
        }
    }
    
    file << "ENDMDL\n";
    file.close();
}

// Save trajectory
void GCMCModule::saveTrajectory(const std::string& filename) {
    config_.trajectoryFile = filename;
    saveSnapshot();
}

// Calculate system energy
double GCMCModule::calculateSystemEnergy() {
    return engine_->calculateSystemEnergy();
}

// Calculate fragment energy
double GCMCModule::calculateFragmentEnergy(int residueIdx) {
    return engine_->calculateFragmentEnergy(residueIdx);
}

// Calculate energy components
std::pair<double, double> GCMCModule::calculateEnergyComponents() {
    // Sum up VDW and electrostatic components
    double vdwEnergy = 0.0;
    double elecEnergy = 0.0;
    
    for (const auto& residue : state_->residues) {
        if (residue.active) {
            vdwEnergy += residue.energy_vdw;
            elecEnergy += residue.energy_elec;
        }
    }
    
    return {elecEnergy, vdwEnergy};
}

// Update statistics
void GCMCModule::updateStatistics() {
    // Get current molecule counts
    std::vector<int> moleculeCounts;
    for (int typeId = 0; typeId < reservoir_->getTemplateCount(); ++typeId) {
        moleculeCounts.push_back(reservoir_->getActiveCount(typeId));
    }
    
    // Calculate current energy
    double totalEnergy = calculateSystemEnergy();
    auto [elecEnergy, vdwEnergy] = calculateEnergyComponents();
    
    // Record state
    statistics_->recordState(moleculeCounts, totalEnergy, vdwEnergy, elecEnergy);
}

// Check convergence
void GCMCModule::checkConvergence() {
    // Check if all fragment types have converged
    bool converged = true;
    for (int typeId = 0; typeId < reservoir_->getTemplateCount(); ++typeId) {
        if (!statistics_->isNumberConverged(typeId, 0.01)) {
            converged = false;
            break;
        }
    }
    
    // Check energy convergence
    if (!statistics_->isEnergyConverged(0.01)) {
        converged = false;
    }
    
    if (converged && config_.verbose) {
        std::cout << "System has converged!" << std::endl;
    }
}

// Adapt biasing
void GCMCModule::adaptBiasing() {
    if (moveSelector_) {
        moveSelector_->updateAdaptiveProbabilities();
    }
    
    if (biasCalc_) {
        biasCalc_->updateBiasParameters();
    }
}

// Enable adaptive biasing
void GCMCModule::enableAdaptiveBiasing() {
    if (moveSelector_) {
        moveSelector_->enableAdaptiveProbabilities();
    }
    
    if (biasCalc_) {
        biasCalc_->enableAdaptiveBiasing();
    }
}

// Set target density
void GCMCModule::setTargetDensity(double density) {
    // Calculate target number of molecules from density
    if (state_ && state_->periodicBox.size() == 3) {
        double volume = state_->periodicBox[0] * state_->periodicBox[1] * state_->periodicBox[2];
        double targetN = density * volume * 6.022e23 / 1e24;  // Convert to molecules/nm^3
        
        // Adjust chemical potentials to achieve target density
        // This is a placeholder - actual implementation would use iterative adjustment
        if (config_.verbose) {
            std::cout << "Target density set to " << density << " g/cm^3" << std::endl;
            std::cout << "Target N = " << targetN << " molecules" << std::endl;
        }
    }
}

// Enable flat histogram
void GCMCModule::enableFlatHistogram() {
    // Placeholder for Wang-Landau or similar methods
    if (config_.verbose) {
        std::cout << "Flat histogram sampling enabled (not yet implemented)" << std::endl;
    }
}

// Validate state consistency (debug only)
void GCMCModule::validateStateConsistency() {
    if (!state_ || !reservoir_) return;
    
    // Check that all active instances in reservoir have corresponding active residues
    auto activeInstances = reservoir_->getActiveInstances();
    for (int instanceId : activeInstances) {
        // Check bounds
        if (instanceId < 0 || instanceId >= static_cast<int>(state_->residues.size())) {
            std::cerr << "ERROR: Active instance " << instanceId 
                     << " out of bounds (residues.size=" << state_->residues.size() << ")" << std::endl;
            assert(false);
        }
        
        // Check if residue is active
        if (!state_->residues[instanceId].active) {
            std::cerr << "ERROR: Instance " << instanceId 
                     << " is active in reservoir but inactive in MCState" << std::endl;
            assert(false);
        }
        
        // CRITICAL: Verify instance-residue index consistency
        FragmentInstance* instance = reservoir_->getInstance(instanceId);
        if (instance) {
            if (instance->instanceId != instanceId) {
                std::cerr << "ERROR: Instance ID mismatch: expected " << instanceId 
                         << " got " << instance->instanceId << std::endl;
                assert(false);
            }
            if (instance->residueIndex != instanceId) {
                std::cerr << "ERROR: Instance " << instanceId 
                         << " has residueIndex=" << instance->residueIndex 
                         << " (should match instanceId)" << std::endl;
                assert(false);
            }
            
            // Verify atoms exist for active residue
            const FragmentTemplate* tmpl = reservoir_->getTemplate(instance->templateId);
            if (tmpl && state_->residues[instanceId].atoms.size() != tmpl->atoms.size()) {
                std::cerr << "ERROR: Residue " << instanceId 
                         << " has " << state_->residues[instanceId].atoms.size() 
                         << " atoms but template has " << tmpl->atoms.size() << std::endl;
                assert(false);
            }
        }
    }
    
    // Check that inactive residues don't have atoms
    for (size_t i = 0; i < state_->residues.size(); ++i) {
        if (!state_->residues[i].active && !state_->residues[i].atoms.empty()) {
            std::cerr << "ERROR: Inactive residue " << i << " has " 
                     << state_->residues[i].atoms.size() << " atoms" << std::endl;
            assert(false);
        }
    }
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc