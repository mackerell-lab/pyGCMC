#include "GCMCSimulation.hpp"
#include "../../movement/reservoir/MultiTypeReservoir.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

using namespace movement::gcmc;

// Constants
constexpr double kB = 8.314462618e-3;  // Boltzmann constant in kJ/(mol*K)
constexpr double NA = 6.02214076e23;   // Avogadro's number

GCMCSimulation::GCMCSimulation(const Config& config) 
    : config_(config), uniform_(0.0, 1.0) {
    
    // Initialize random number generator
    if (config.randomSeed < 0) {
        std::random_device rd;
        rng_.seed(rd());
    } else {
        rng_.seed(config.randomSeed);
    }
    
    // Set up logging
    if (config.verbose) {
        system::log::LogMain::set_verbose(true);
        // Set verbose logging
        // system::log::LogMain::set_log_level(system::log::LogLevel::DEBUG);
    }
}

GCMCSimulation::~GCMCSimulation() {
    if (initialized_) {
        finalize();
    }
}

bool GCMCSimulation::initialize() {
    log("Initializing GCMC simulation from %s", config_.inputFile.c_str());
    
    // Load parameters from input file
    if (!loadParameters()) {
        log("ERROR: Failed to load parameters");
        return false;
    }
    // If print frequency wasn't explicitly set via CLI, use INP nprint
    if (config_.printFrequency <= 0) {
        config_.printFrequency = params_->get_mc_info().print_freq;
    }
    
    // Setup the MC state
    if (!setupSystem()) {
        log("ERROR: Failed to setup system");
        return false;
    }
    
    // Setup fragment types and reservoir
    if (!setupFragments()) {
        log("ERROR: Failed to setup fragments");
        return false;
    }
    
    // Setup acceptance calculator
    if (!setupAcceptance()) {
        log("ERROR: Failed to setup acceptance calculator");
        return false;
    }
    
    // Setup GCMC engine
    if (!setupEngine()) {
        log("ERROR: Failed to setup GCMC engine");
        return false;
    }
    
    // Initialize statistics
    statistics_.setAutoAdjust(config_.enableAdaptiveSampling);
    statistics_.setSamplingInterval(config_.statisticsInterval);
    
    initialized_ = true;
    log("Initialization complete");
    
    // Print initial system information
    printStatistics();
    
    return true;
}

bool GCMCSimulation::loadParameters() {
    try {
        params_ = std::make_unique<model::param::Param>();
        
        // Parse INP file (use extended GCMC parser)
        io::parameters::InpParserGCMC::parse_to_param(config_.inputFile, *params_);
        
        // Update derived values
        params_->update_derived_values();
        
        log("Loaded parameters from %s", config_.inputFile.c_str());
        log("  Temperature: %.2f K", params_->get_mc_info().temperature);
        log("  Box size: %.2f x %.2f x %.2f nm", 
            params_->get_space_info().box_size[0],
            params_->get_space_info().box_size[1],
            params_->get_space_info().box_size[2]);
        log("  MC steps: %d", params_->get_mc_info().mc_steps);
        
        return true;
        
    } catch (const std::exception& e) {
        log("ERROR: Exception loading parameters: %s", e.what());
        return false;
    }
}

bool GCMCSimulation::setupSystem() {
    state_ = std::make_unique<model::montecarlo::MCState>();
    
    // Set box dimensions
    const auto& box = params_->get_space_info().box_size;
    state_->info.box[0] = box[0];
    state_->info.box[1] = box[1];
    state_->info.box[2] = box[2];
    
    // Set temperature through beta from parameters
    const double beta = params_->get_mc_info().beta;
    state_->info.beta = beta;
    
    // Load initial structure if provided
    const auto& pdbFile = params_->get_file_info().input_pdb_file;
    if (!pdbFile.empty()) {
        log("Loading initial structure from %s", pdbFile.c_str());
        // TODO: Integrate PDB parser when available
        /*
        try {
            auto structure = io::PdbParserMain::parse_file(pdbFile);
            
            // Convert PDB structure to MCState atoms and residues
            state_->atoms.resize(structure.get_atoms().size());
            
            for (size_t i = 0; i < structure.get_atoms().size(); ++i) {
                const auto& pdbAtom = structure.get_atoms()[i];
                auto& mcAtom = state_->atoms[i];
                mcAtom.x = pdbAtom.get_x() / 10.0;  // Convert Angstrom to nm
                mcAtom.y = pdbAtom.get_y() / 10.0;
                mcAtom.z = pdbAtom.get_z() / 10.0;
                mcAtom.type = 0;  // Will be assigned from topology
                mcAtom.name = pdbAtom.get_name();
                mcAtom.updatePosition();
            }
            
            // Convert residues
            state_->residues.resize(structure.get_residues().size());
            
            for (size_t i = 0; i < structure.get_residues().size(); ++i) {
                const auto& pdbRes = structure.get_residues()[i];
                auto& mcRes = state_->residues[i];
                mcRes.resname = pdbRes.get_name();
                mcRes.resid = pdbRes.get_resseq();
                mcRes.atomStart = pdbRes.get_atom_indices().empty() ? 0 : pdbRes.get_atom_indices()[0];
                mcRes.atomCount = pdbRes.get_atom_indices().size();
                mcRes.active = true;
                mcRes.fixed = false;
                
                // Calculate center of mass
                mcRes.center[0] = pdbRes.get_center_of_mass()[0] / 10.0;  // Convert to nm
                mcRes.center[1] = pdbRes.get_center_of_mass()[1] / 10.0;
                mcRes.center[2] = pdbRes.get_center_of_mass()[2] / 10.0;
            }
            
            log("Loaded %zu atoms and %zu residues from PDB", 
                state_->atoms.size(), state_->residues.size());
                
        } catch (const std::exception& e) {
            log("ERROR: Failed to parse PDB file: %s", e.what());
            return false;
        }
        */
        log("PDB loading temporarily disabled - using empty initial state");
    }
    
    // Load topology and setup force field
    const auto& topFile = params_->get_file_info().topology_file;
    if (!topFile.empty()) {
        log("Loading topology from %s", topFile.c_str());
        // TODO: Integrate TOP parser when available
        /*
        try {
            auto topology = io::TOPParser::parse_file(topFile);
            
            // Extract atom types and build force field
            // For now, use simple placeholder values
            // TODO: Integrate with proper force field parameters from PAR files
            size_t numTypes = 10;  // Placeholder - should come from topology
            state_->forcefield.numTotalTypes = numTypes;
            state_->forcefield.numMovementTypes = 4;  // Water, ions, etc.
            
            // Initialize LJ parameters with placeholder values
            state_->forcefield.ljSigma.resize(numTypes * numTypes);
            state_->forcefield.ljEpsilon.resize(numTypes * numTypes);
            
            for (size_t i = 0; i < numTypes; ++i) {
                for (size_t j = 0; j < numTypes; ++j) {
                    size_t idx = i * numTypes + j;
                    // Placeholder LJ parameters (will be replaced with real values from PAR files)
                    state_->forcefield.ljSigma[idx] = 0.3f + 0.01f * (i + j);  // nm
                    state_->forcefield.ljEps[idx] = 0.5f + 0.05f * (i * j);  // kJ/mol
                }
            }
            
            log("Initialized force field with %zu types", numTypes);
            
        } catch (const std::exception& e) {
            log("WARNING: Failed to parse topology file: %s", e.what());
            log("Using default force field parameters");
            // Fall back to placeholder values
            state_->forcefield.numTotalTypes = 10;
            state_->forcefield.numMovementTypes = 4;
        }
        */
        log("Topology loading temporarily disabled - using placeholder force field");
    }
    
    // Use default force field for now
    size_t numTypes = 10;  // Placeholder
    state_->forcefield.numTotalTypes = numTypes;
    state_->forcefield.numMovementTypes = 4;
    
    // Initialize LJ parameters with placeholder values
    state_->forcefield.ljSigma.resize(numTypes * numTypes);
    state_->forcefield.ljEps.resize(numTypes * numTypes);
    
    for (size_t i = 0; i < numTypes; ++i) {
        for (size_t j = 0; j < numTypes; ++j) {
            size_t idx = i * numTypes + j;
            state_->forcefield.ljSigma[idx] = 0.3f + 0.01f * (i + j);  // nm
            state_->forcefield.ljEps[idx] = 0.5f + 0.05f * (i * j);  // kJ/mol
        }
    }
    
    return true;
}

bool GCMCSimulation::setupFragments() {
    const auto& fragInfo = params_->get_fragment_info();
    const auto& fileInfo = params_->get_file_info();
    const auto& mcInfo = params_->get_mc_info();
    
    // Create multi-type fragment reservoir
    reservoir_ = std::make_unique<movement::MultiTypeReservoir>();
    
    // Process each fragment type
    for (size_t i = 0; i < fileInfo.fragment_names.size(); ++i) {
        FragmentInfo frag;
        frag.name = fileInfo.fragment_names[i];
        frag.typeId = i;
        
        // Set concentration and chemical potential
        if (i < fragInfo.conc_list.size()) {
            frag.concentration = fragInfo.conc_list[i];
        }
        if (i < fragInfo.muex_list.size()) {
            frag.chemicalPotential = fragInfo.muex_list[i];
        }
        
        // Calculate activity from chemical potential
        double beta = state_->info.beta;
        frag.activity = std::exp(beta * frag.chemicalPotential);
        
        // Set probability from MC time allocation
        if (i < mcInfo.mc_time_list.size()) {
            frag.probability = mcInfo.mc_time_list[i];
        } else {
            frag.probability = 1.0 / fileInfo.fragment_names.size();
        }
        
        // Calculate maximum count from concentration and volume
        double volume = params_->get_space_info().volume;
        // Convert concentration (M) to number: N = C * V * NA / 1000
        frag.maxCount = static_cast<int>(frag.concentration * volume * NA / 1000.0);
        
        // Create fragment template
        movement::FragmentTemplate tmpl;
        tmpl.name = frag.name;
        tmpl.typeId = frag.typeId;
        tmpl.chemicalPotential = frag.chemicalPotential;
        tmpl.activity = frag.activity;
        tmpl.concentration = frag.concentration;
        tmpl.radius = (i < fragInfo.radius_list.size()) ? fragInfo.radius_list[i] : 0.0;
        
        // Create fragment atoms based on fragment type
        // TODO: Replace with ITP parser when ready
        if (frag.name == "water" || frag.name == "WAT" || frag.name == "HOH") {
            // Water molecule: O-H-H
            tmpl.atoms.resize(3);
            // Oxygen
            tmpl.atoms[0].x = 0.0;
            tmpl.atoms[0].y = 0.0;
            tmpl.atoms[0].z = 0.0;
            tmpl.atoms[0].type = 0;  // O type
            tmpl.atoms[0].charge = -0.834;  // TIP3P charge
            tmpl.atoms[0].mass = 15.999;
            // H1
            tmpl.atoms[1].x = 0.0756;
            tmpl.atoms[1].y = 0.0586;
            tmpl.atoms[1].z = 0.0;
            tmpl.atoms[1].type = 1;  // H type
            tmpl.atoms[1].charge = 0.417;
            tmpl.atoms[1].mass = 1.008;
            // H2
            tmpl.atoms[2].x = -0.0756;
            tmpl.atoms[2].y = 0.0586;
            tmpl.atoms[2].z = 0.0;
            tmpl.atoms[2].type = 1;  // H type
            tmpl.atoms[2].charge = 0.417;
            tmpl.atoms[2].mass = 1.008;
        } else if (frag.name == "Na" || frag.name == "NA" || frag.name == "SOD") {
            // Sodium ion
            tmpl.atoms.resize(1);
            tmpl.atoms[0].x = 0.0;
            tmpl.atoms[0].y = 0.0;
            tmpl.atoms[0].z = 0.0;
            tmpl.atoms[0].type = 2;  // Na type
            tmpl.atoms[0].charge = 1.0;
            tmpl.atoms[0].mass = 22.990;
        } else if (frag.name == "Cl" || frag.name == "CL" || frag.name == "CLA") {
            // Chloride ion
            tmpl.atoms.resize(1);
            tmpl.atoms[0].x = 0.0;
            tmpl.atoms[0].y = 0.0;
            tmpl.atoms[0].z = 0.0;
            tmpl.atoms[0].type = 3;  // Cl type
            tmpl.atoms[0].charge = -1.0;
            tmpl.atoms[0].mass = 35.453;
        } else {
            // Default: single atom placeholder
            tmpl.atoms.resize(1);
            tmpl.atoms[0].x = 0.0;
            tmpl.atoms[0].y = 0.0;
            tmpl.atoms[0].z = 0.0;
            tmpl.atoms[0].type = i % state_->forcefield.numMovementTypes;
            tmpl.atoms[0].charge = 0.0;
            tmpl.atoms[0].mass = 12.0;
        }
        
        // Create type info for multi-type reservoir
        movement::MultiTypeReservoir::TypeInfo typeInfo;
        typeInfo.typeId = frag.typeId;
        typeInfo.name = frag.name;
        typeInfo.chemicalPotential = frag.chemicalPotential;
        typeInfo.activity = frag.activity;
        typeInfo.probability = frag.probability;
        typeInfo.maxCount = frag.maxCount;
        typeInfo.radius = tmpl.radius;
        
        // Add to reservoir
        reservoir_->addType(typeInfo, tmpl);
        
        // Store fragment info
        fragmentTypes_.push_back(frag);
        fragmentNameToId_[frag.name] = frag.typeId;
        
        log("Fragment %s: conc=%.3f M, mu=%.2f kJ/mol, activity=%.3e, max=%d",
            frag.name.c_str(), frag.concentration, frag.chemicalPotential,
            frag.activity, frag.maxCount);
    }
    
    return true;
}

bool GCMCSimulation::setupAcceptance() {
    acceptance_ = std::make_unique<GCMCAcceptance>();
    
    // Set temperature from parameters
    const double temperature = params_->get_mc_info().temperature;
    acceptance_->setTemperature(temperature);
    
    // Set volume
    const auto& box = state_->info.box;
    double volume = box[0] * box[1] * box[2];
    acceptance_->setVolume(volume);
    
    // Set activities for each fragment type
    for (const auto& frag : fragmentTypes_) {
        acceptance_->setActivity(frag.typeId, frag.activity);
    }
    
    log("Acceptance calculator configured:");
    log("  Temperature: %.2f K", params_->get_mc_info().temperature);
    log("  Volume: %.2f nm^3", volume);
    
    return true;
}

bool GCMCSimulation::setupEngine() {
    engine_ = std::make_unique<GCMCEngine>();
    
    // Initialize with state and reservoir
    engine_->initialize(state_.get(), reservoir_.get());
    
    // Set acceptance calculator
    engine_->setAcceptanceCalculator(acceptance_.get());
    
    // Configure engine parameters
    engine_->setConfigValue("maxTranslation", 0.2);  // nm
    engine_->setConfigValue("maxRotation", 15.0);    // degrees
    engine_->setConfigValue("useCavityBias", params_->get_bias_info().use_cavity_bias ? 1.0 : 0.0);
    
    // Enable statistics if configured
    if (config_.enableStatistics) {
        engine_->enableStatistics(true);
        engine_->setStatisticsInterval(config_.statisticsInterval);
    }
    
    // Set probability storage
    if (config_.storeProbabilities) {
        engine_->setConfigValue("storeProbabilities", 1.0);
    }
    
    log("GCMC engine configured:");
    log("  Max translation: %.3f nm", engine_->getConfigValue("maxTranslation"));
    log("  Max rotation: %.1f degrees", engine_->getConfigValue("maxRotation"));
    log("  Cavity bias: %s", params_->get_bias_info().use_cavity_bias ? "enabled" : "disabled");
    
    // Ensure deterministic RNG when seed provided
    if (config_.randomSeed >= 0) {
        engine_->setSeed(static_cast<unsigned int>(config_.randomSeed));
    }

    return true;
}

bool GCMCSimulation::run() {
    if (!initialized_) {
        log("ERROR: Simulation not initialized");
        return false;
    }
    
    log("Starting GCMC simulation for %d steps", params_->get_mc_info().mc_steps);
    
    running_ = true;
    startTime_ = std::chrono::steady_clock::now();
    
    int mcSteps = params_->get_mc_info().mc_steps;
    
    for (int step = 0; step < mcSteps && running_; ++step) {
        // Perform MC move
        if (!performMCStep()) {
            log("ERROR: Failed at step %d", step);
            return false;
        }
        
        stats_.totalSteps++;
        
        // Print statistics
        if (step % config_.printFrequency == 0 && step > 0) {
            writeStatistics(step);
        }
        
        // Save trajectory
        if (step % config_.trajectoryFrequency == 0 && step > 0) {
            writeTrajectory(step);
        }
        
        // Save checkpoint
        if (step % config_.checkpointFrequency == 0 && step > 0) {
            writeCheckpoint(step);
        }
        
        // Check convergence
        if (config_.enableAdaptiveSampling && step % 10000 == 0 && step > 0) {
            if (checkConvergence()) {
                log("Simulation converged at step %d", step);
                break;
            }
        }
    }
    
    running_ = false;
    
    // Calculate final statistics
    auto endTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime_;
    stats_.totalTime = elapsed.count();
    stats_.timePerStep = stats_.totalTime / stats_.totalSteps;
    stats_.stepsPerSecond = stats_.totalSteps / stats_.totalTime;
    
    log("Simulation completed:");
    log("  Total steps: %d", stats_.totalSteps);
    log("  Total time: %.2f seconds", stats_.totalTime);
    log("  Performance: %.1f steps/second", stats_.stepsPerSecond);
    // Emit a completion line to stdout even in non-verbose mode (tests expect this)
    std::cout << "Simulation completed" << std::endl;
    
    return true;
}

bool GCMCSimulation::performMCStep() {
    // Select move type
    MoveType moveType = selectMoveType();
    
    bool accepted = false;
    GCMCEngine::MoveResult result;
    
    switch (moveType) {
        case INSERT: {
            int fragType = selectFragmentType();
            result = engine_->attemptInsertion(fragType);
            
            fragmentTypes_[fragType].insertAttempts++;
            if (result.accepted) {
                fragmentTypes_[fragType].insertAccepted++;
                fragmentTypes_[fragType].currentCount++;
                accepted = true;
            }
            
            stats_.moveAttempts["insertion"]++;
            if (accepted) stats_.moveAccepted["insertion"]++;
            
            // Update new statistics module
            simulationStats_.recordMove("insert", fragmentTypes_[fragType].name, accepted);
            break;
        }
        
        case DELETE: {
            int fragType = -1;
            if (reservoir_->getActiveCount() > 0) {
                fragType = selectActiveFragment();
                if (fragType >= 0) {
                    result = engine_->attemptDeletion(fragType);
                    
                    fragmentTypes_[fragType].deleteAttempts++;
                    if (result.accepted) {
                        fragmentTypes_[fragType].deleteAccepted++;
                        fragmentTypes_[fragType].currentCount--;
                        accepted = true;
                    }
                }
            }
            
            stats_.moveAttempts["deletion"]++;
            if (accepted) stats_.moveAccepted["deletion"]++;
            
            // Update new statistics module
            if (fragType >= 0) {
                simulationStats_.recordMove("delete", fragmentTypes_[fragType].name, accepted);
            }
            break;
        }
        
        case TRANSLATE: {
            std::vector<int> activeIndices;
            if (reservoir_->getActiveCount() > 0) {
                activeIndices = reservoir_->getActiveInstances();
                if (!activeIndices.empty()) {
                    int idx = activeIndices[rng_() % activeIndices.size()];
                    result = engine_->attemptTranslation(idx);
                    accepted = result.accepted;
                    
                    // Update new statistics module for translation
                    // Note: We'll skip type lookup for now as getInstanceType doesn't exist
                    // Just record with generic "fragment" name
                    simulationStats_.recordMove("translate", "fragment", accepted);
                }
            }
            
            stats_.moveAttempts["translation"]++;
            if (accepted) stats_.moveAccepted["translation"]++;
            break;
        }
        
        case ROTATE: {
            std::vector<int> activeIndices;
            if (reservoir_->getActiveCount() > 0) {
                activeIndices = reservoir_->getActiveInstances();
                if (!activeIndices.empty()) {
                    int idx = activeIndices[rng_() % activeIndices.size()];
                    result = engine_->attemptRotation(idx);
                    accepted = result.accepted;
                    
                    // Update new statistics module for rotation
                    // Note: We'll skip type lookup for now as getInstanceType doesn't exist
                    // Just record with generic "fragment" name
                    simulationStats_.recordMove("rotate", "fragment", accepted);
                }
            }
            
            stats_.moveAttempts["rotation"]++;
            if (accepted) stats_.moveAccepted["rotation"]++;
            break;
        }
    }
    
    if (accepted) {
        stats_.acceptedMoves++;
    }
    
    // Update energy history
    if (config_.enableStatistics && stats_.totalSteps % config_.statisticsInterval == 0) {
        double energy = calculateSystemEnergy();
        stats_.energyHistory.push_back(energy);
        stats_.currentEnergy = energy;
        
        // Record in new statistics module
        simulationStats_.recordEnergy(energy);
    }
    
    return true;
}

GCMCSimulation::MoveType GCMCSimulation::selectMoveType() {
    double r = uniform_(rng_);
    
    // Simple equal probability for now
    // TODO: Implement adaptive move probabilities
    if (r < 0.25) return INSERT;
    else if (r < 0.50) return DELETE;
    else if (r < 0.75) return TRANSLATE;
    else return ROTATE;
}

int GCMCSimulation::selectFragmentType() {
    if (fragmentTypes_.size() == 1) {
        return 0;
    }
    
    // Use weighted selection based on MC time allocation
    double r = uniform_(rng_);
    double cumSum = 0.0;
    
    for (size_t i = 0; i < fragmentTypes_.size(); ++i) {
        cumSum += fragmentTypes_[i].probability;
        if (r < cumSum) {
            return i;
        }
    }
    
    return fragmentTypes_.size() - 1;
}

int GCMCSimulation::selectActiveFragment() {
    // For deletion, select uniformly from active fragments
    auto activeIndices = reservoir_->getActiveInstances();
    if (activeIndices.empty()) {
        return -1;
    }
    
    int idx = activeIndices[rng_() % activeIndices.size()];
    auto instance = reservoir_->getInstance(idx);
    if (instance && instance->isActive) {
        // Get type from template
        const auto* tpl = reservoir_->getTemplate(instance->templateId);
        return tpl ? tpl->typeId : -1;
    }
    
    return -1;
}

double GCMCSimulation::calculateSystemEnergy() {
    // Use the engine's energy calculation
    return engine_->calculateSystemEnergy();
}

void GCMCSimulation::updateStatistics() {
    // Update acceptance rates
    stats_.acceptanceRate = static_cast<double>(stats_.acceptedMoves) / stats_.totalSteps;
    
    for (auto& [move, attempts] : stats_.moveAttempts) {
        if (attempts > 0) {
            stats_.moveAcceptanceRates[move] = 
                static_cast<double>(stats_.moveAccepted[move]) / attempts;
        }
    }
    
    // Update fragment statistics
    for (auto& frag : fragmentTypes_) {
        stats_.fragmentCounts[frag.name] = frag.currentCount;
        
        // Calculate density (molecules/nm^3)
        double volume = state_->info.box[0] * state_->info.box[1] * state_->info.box[2];
        stats_.fragmentDensities[frag.name] = frag.currentCount / volume;
        
        // Calculate acceptance rates
        if (frag.insertAttempts > 0) {
            double insertRate = static_cast<double>(frag.insertAccepted) / frag.insertAttempts;
            double deleteRate = (frag.deleteAttempts > 0) ? 
                static_cast<double>(frag.deleteAccepted) / frag.deleteAttempts : 0.0;
            stats_.fragmentAcceptanceRates[frag.name] = (insertRate + deleteRate) / 2.0;
        }
    }
    
    // Update energy statistics
    if (!stats_.energyHistory.empty()) {
        double sum = 0.0;
        for (double e : stats_.energyHistory) {
            sum += e;
        }
        stats_.averageEnergy = sum / stats_.energyHistory.size();
        
        // Calculate standard deviation
        double sumSq = 0.0;
        for (double e : stats_.energyHistory) {
            double diff = e - stats_.averageEnergy;
            sumSq += diff * diff;
        }
        stats_.energyStdDev = std::sqrt(sumSq / stats_.energyHistory.size());
    }
}

bool GCMCSimulation::checkConvergence() {
    // Simple convergence check based on energy fluctuations
    if (stats_.energyHistory.size() < 100) {
        return false;
    }
    
    // Check if energy standard deviation is small relative to average
    if (stats_.averageEnergy != 0.0) {
        double relStdDev = stats_.energyStdDev / std::abs(stats_.averageEnergy);
        return relStdDev < config_.convergenceTolerance;
    }
    
    return false;
}

void GCMCSimulation::writeStatistics(int step) {
    updateStatistics();
    
    std::cout << "\n=== Step " << step << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    
    std::cout << "Acceptance: " << stats_.acceptanceRate * 100 << "%" << std::endl;
    std::cout << "Energy: " << stats_.currentEnergy << " kJ/mol" << std::endl;
    
    std::cout << "Fragment counts:" << std::endl;
    for (const auto& frag : fragmentTypes_) {
        std::cout << "  " << frag.name << ": " << frag.currentCount;
        if (frag.insertAttempts > 0) {
            double rate = static_cast<double>(frag.insertAccepted) / frag.insertAttempts;
            std::cout << " (accept: " << rate * 100 << "%)";
        }
        std::cout << std::endl;
    }
    
    std::cout << "Performance: " << stats_.totalSteps / stats_.totalTime << " steps/s" << std::endl;
}

void GCMCSimulation::writeTrajectory(int step) {
    std::string filename = config_.outputPrefix + "_traj_" + std::to_string(step) + ".pdb";
    saveTrajectory(filename);
    log("Saved trajectory to %s", filename.c_str());
}

void GCMCSimulation::writeCheckpoint(int step) {
    std::string filename = config_.outputPrefix + "_checkpoint_" + std::to_string(step) + ".dat";
    saveCheckpoint(filename);
    log("Saved checkpoint to %s", filename.c_str());
}

void GCMCSimulation::writeFinalResults() {
    updateStatistics();
    
    std::string filename = config_.outputPrefix + "_final.txt";
    std::ofstream out(filename);
    
    out << "GCMC Simulation Final Results\n";
    out << "==============================\n\n";
    
    out << "Configuration:\n";
    out << "  Input file: " << config_.inputFile << "\n";
    out << "  Temperature: " << params_->get_mc_info().temperature << " K\n";
    out << "  Box: " << state_->info.box[0] << " x " << state_->info.box[1] 
        << " x " << state_->info.box[2] << " nm\n";
    out << "  Total steps: " << stats_.totalSteps << "\n\n";
    
    if (config_.enableStatistics) {
        out << "Performance:\n";
        out << "  Total time: " << stats_.totalTime << " seconds\n";
        out << "  Steps/second: " << stats_.stepsPerSecond << "\n\n";
    }
    
    out << "Statistics:\n";
    out << "  Overall acceptance: " << stats_.acceptanceRate * 100 << "%\n";
    out << "  Average energy: " << stats_.averageEnergy << " +/- " 
        << stats_.energyStdDev << " kJ/mol\n\n";
    
    out << "Fragment Statistics:\n";
    for (const auto& frag : fragmentTypes_) {
        out << "  " << frag.name << ":\n";
        out << "    Final count: " << frag.currentCount << "\n";
        out << "    Density: " << stats_.fragmentDensities[frag.name] << " molecules/nm^3\n";
        out << "    Insert attempts: " << frag.insertAttempts << "\n";
        out << "    Insert accepted: " << frag.insertAccepted << "\n";
        out << "    Delete attempts: " << frag.deleteAttempts << "\n";
        out << "    Delete accepted: " << frag.deleteAccepted << "\n";
    }
    
    out.close();
    log("Wrote final results to %s", filename.c_str());
}

void GCMCSimulation::finalize() {
    if (!initialized_) return;
    
    writeFinalResults();
    
    // Save final trajectory
    std::string trajFile = config_.outputPrefix + "_final.pdb";
    saveTrajectory(trajFile);
    
    // Save final checkpoint
    std::string checkFile = config_.outputPrefix + "_final.checkpoint";
    saveCheckpoint(checkFile);
    
    log("Simulation finalized");
}

void GCMCSimulation::printStatistics() const {
    // Use the new statistics module for formatted output
    simulationStats_.printSummary(stats_.totalSteps);
    
    // Also print legacy statistics if needed
    if (config_.verbose) {
        std::cout << "\n=== GCMC Simulation Statistics ===" << std::endl;
        std::cout << "Total steps: " << stats_.totalSteps << std::endl;
        std::cout << "Accepted moves: " << stats_.acceptedMoves << std::endl;
        if (stats_.totalSteps > 0) {
            std::cout << "Acceptance rate: " 
                      << (100.0 * stats_.acceptedMoves / stats_.totalSteps) << "%" << std::endl;
        }
        
        for (const auto& frag : fragmentTypes_) {
            std::cout << "Fragment " << frag.name << ": " << frag.currentCount << " molecules" << std::endl;
        }
    }
}

void GCMCSimulation::saveTrajectory(const std::string& filename) const {
    if (!state_) return;
    
    std::ofstream out(filename);
    if (!out) {
        log("ERROR: Failed to open trajectory file %s", filename.c_str());
        return;
    }
    
    // Write PDB header
    out << "REMARK GCMC Trajectory\n";
    out << "REMARK Step: " << stats_.totalSteps << "\n";
    out << "REMARK Energy: " << stats_.currentEnergy << " kJ/mol\n";
    
    // Write box dimensions (CRYST1 record)
    out << "CRYST1";
    out << std::fixed << std::setprecision(3);
    out << std::setw(9) << state_->info.box[0] * 10.0;  // nm to Angstrom
    out << std::setw(9) << state_->info.box[1] * 10.0;
    out << std::setw(9) << state_->info.box[2] * 10.0;
    out << std::setw(7) << "90.00";
    out << std::setw(7) << "90.00";
    out << std::setw(7) << "90.00";
    out << " P 1           1\n";
    
    // Write atoms
    int atomIdx = 1;
    for (size_t resIdx = 0; resIdx < state_->residues.size(); ++resIdx) {
        const auto& res = state_->residues[resIdx];
        if (!res.active) continue;
        
        // Write atoms for this residue
        for (int j = 0; j < res.atomCount && (res.atomStart + j) < static_cast<int>(state_->atoms.size()); ++j) {
            const auto& atom = state_->atoms[res.atomStart + j];
            
            out << "ATOM  ";
            out << std::setw(5) << atomIdx++;
            out << "  ";
            
            // Atom name - based on type
            if (atom.type == 0) {
                out << std::left << std::setw(4) << "O";
            } else if (atom.type == 1) {
                out << std::left << std::setw(4) << "H";
            } else if (atom.type == 2) {
                out << std::left << std::setw(4) << "Na";
            } else if (atom.type == 3) {
                out << std::left << std::setw(4) << "Cl";
            } else {
                out << std::left << std::setw(4) << "X";
            }
            
            // Residue name and number
            out << std::right;
            out << std::setw(3) << res.resname.substr(0, 3);
            out << " A";  // Default chain ID
            out << std::setw(4) << res.resid;
            out << "    ";
            
            // Coordinates (nm to Angstrom)
            out << std::fixed << std::setprecision(3);
            out << std::setw(8) << atom.x * 10.0;
            out << std::setw(8) << atom.y * 10.0;
            out << std::setw(8) << atom.z * 10.0;
            
            // Occupancy and temperature factor
            out << std::setw(6) << "1.00";
            out << std::setw(6) << "0.00";
            
            // Element symbol
            out << "          ";
            if (atom.type == 0) out << " O";
            else if (atom.type == 1) out << " H";
            else if (atom.type == 2) out << "Na";
            else if (atom.type == 3) out << "Cl";
            else out << " X";
            
            out << "\n";
        }
    }
    
    out << "END\n";
    out.close();
    
    log("Saved trajectory to %s", filename.c_str());
}

void GCMCSimulation::saveCheckpoint(const std::string& filename) const {
    if (!state_ || !engine_) return;
    
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        log("ERROR: Failed to open checkpoint file %s", filename.c_str());
        return;
    }
    
    // Write checkpoint header
    const std::string header = "GCMC_CHECKPOINT_V1";
    out.write(header.c_str(), header.size());
    
    // Write simulation state
    out.write(reinterpret_cast<const char*>(&stats_.totalSteps), sizeof(stats_.totalSteps));
    out.write(reinterpret_cast<const char*>(&stats_.acceptedMoves), sizeof(stats_.acceptedMoves));
    out.write(reinterpret_cast<const char*>(&stats_.currentEnergy), sizeof(stats_.currentEnergy));
    out.write(reinterpret_cast<const char*>(&stats_.totalTime), sizeof(stats_.totalTime));
    
    // Write fragment counts
    size_t numFragTypes = fragmentTypes_.size();
    out.write(reinterpret_cast<const char*>(&numFragTypes), sizeof(numFragTypes));
    for (const auto& frag : fragmentTypes_) {
        size_t nameLen = frag.name.size();
        out.write(reinterpret_cast<const char*>(&nameLen), sizeof(nameLen));
        out.write(frag.name.c_str(), nameLen);
        out.write(reinterpret_cast<const char*>(&frag.currentCount), sizeof(frag.currentCount));
        out.write(reinterpret_cast<const char*>(&frag.insertAttempts), sizeof(frag.insertAttempts));
        out.write(reinterpret_cast<const char*>(&frag.insertAccepted), sizeof(frag.insertAccepted));
        out.write(reinterpret_cast<const char*>(&frag.deleteAttempts), sizeof(frag.deleteAttempts));
        out.write(reinterpret_cast<const char*>(&frag.deleteAccepted), sizeof(frag.deleteAccepted));
    }
    
    // Write atom positions
    size_t numAtoms = state_->atoms.size();
    out.write(reinterpret_cast<const char*>(&numAtoms), sizeof(numAtoms));
    for (const auto& atom : state_->atoms) {
        out.write(reinterpret_cast<const char*>(&atom.x), sizeof(atom.x));
        out.write(reinterpret_cast<const char*>(&atom.y), sizeof(atom.y));
        out.write(reinterpret_cast<const char*>(&atom.z), sizeof(atom.z));
        out.write(reinterpret_cast<const char*>(&atom.type), sizeof(atom.type));
        // MCAtom doesn't have isActive, write a placeholder
        bool active = true;
        out.write(reinterpret_cast<const char*>(&active), sizeof(active));
    }
    
    // Write residue information
    size_t numResidues = state_->residues.size();
    out.write(reinterpret_cast<const char*>(&numResidues), sizeof(numResidues));
    for (const auto& res : state_->residues) {
        size_t nameLen = res.resname.size();
        out.write(reinterpret_cast<const char*>(&nameLen), sizeof(nameLen));
        out.write(res.resname.c_str(), nameLen);
        // Write active status
        out.write(reinterpret_cast<const char*>(&res.active), sizeof(res.active));
    }
    
    out.close();
    log("Saved checkpoint to %s", filename.c_str());
}

bool GCMCSimulation::loadCheckpoint(const std::string& /*filename*/) {
    // TODO: Implement checkpoint loading
    // Will deserialize state_, reservoir_, and statistics
    return false;
}

template<typename... Args>
void GCMCSimulation::log(const std::string& format, Args... args) const {
    if (config_.verbose) {
        system::log::LogMain::info(format, args...);
    }
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc
