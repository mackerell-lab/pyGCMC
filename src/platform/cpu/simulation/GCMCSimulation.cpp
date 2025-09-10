#include "GCMCSimulation.hpp"
#include "../movement/reservoir/MultiTypeReservoir.hpp"
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
        // TODO: Load PDB file into state
        log("Loading initial structure from %s", pdbFile.c_str());
    }
    
    // Setup force field
    // TODO: Load force field parameters from files
    state_->forcefield.numTotalTypes = 10;  // Placeholder
    state_->forcefield.numMovementTypes = 4;  // Placeholder
    
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
        tmpl.atoms.resize(1);  // Placeholder until ITP parser is ready
        
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
            break;
        }
        
        case DELETE: {
            if (reservoir_->getActiveCount() > 0) {
                int fragType = selectActiveFragment();
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
            break;
        }
        
        case TRANSLATE: {
            if (reservoir_->getActiveCount() > 0) {
                auto activeIndices = reservoir_->getActiveInstances();
                if (!activeIndices.empty()) {
                    int idx = activeIndices[rng_() % activeIndices.size()];
                    result = engine_->attemptTranslation(idx);
                    accepted = result.accepted;
                }
            }
            
            stats_.moveAttempts["translation"]++;
            if (accepted) stats_.moveAccepted["translation"]++;
            break;
        }
        
        case ROTATE: {
            if (reservoir_->getActiveCount() > 0) {
                auto activeIndices = reservoir_->getActiveInstances();
                if (!activeIndices.empty()) {
                    int idx = activeIndices[rng_() % activeIndices.size()];
                    result = engine_->attemptRotation(idx);
                    accepted = result.accepted;
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
    if (stats_.totalSteps % config_.statisticsInterval == 0) {
        double energy = calculateSystemEnergy();
        stats_.energyHistory.push_back(energy);
        stats_.currentEnergy = energy;
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
    
    out << "Performance:\n";
    out << "  Total time: " << stats_.totalTime << " seconds\n";
    out << "  Steps/second: " << stats_.stepsPerSecond << "\n\n";
    
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

void GCMCSimulation::saveTrajectory(const std::string& /*filename*/) const {
    // TODO: Implement PDB trajectory writing
    // Will use state_->residues to write current configuration
}

void GCMCSimulation::saveCheckpoint(const std::string& /*filename*/) const {
    // TODO: Implement checkpoint saving
    // Will serialize state_, reservoir_, and statistics
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