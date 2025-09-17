#include "GCMCSimulationModular.hpp"
#include "../movement/reservoir/fragment_reservoir.hpp"
#include <iostream>
#include <random>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

GCMCSimulationModular::GCMCSimulationModular(const Config& config)
    : config_(config)
    , initialized_(false)
    , running_(false) {
    
    // Create modules
    core_ = std::make_unique<SimulationRunner>(nullptr);  // Will set parent later if needed
    setup_ = std::make_unique<SystemInitializer>();
    
    // Configure IO modules
    io::output::TrajectoryWriter::Config trajConfig;
    trajConfig.prefix = config.outputPrefix;
    trajectoryWriter_ = std::make_unique<io::output::TrajectoryWriter>(trajConfig);
    
    io::output::CheckpointManager::Config chkConfig;
    chkConfig.prefix = config.outputPrefix;
    checkpointManager_ = std::make_unique<io::output::CheckpointManager>(chkConfig);
    
    // Create core components
    params_ = std::make_unique<model::param::Param>();
    state_ = std::make_unique<model::montecarlo::MCState>();
    engine_ = std::make_unique<movement::gcmc::GCMCEngine>();
    acceptance_ = std::make_unique<movement::gcmc::GCMCAcceptance>();
    reservoir_ = std::make_unique<movement::MultiTypeReservoir>();
}

GCMCSimulationModular::~GCMCSimulationModular() {
    if (running_) {
        stop();
    }
}

bool GCMCSimulationModular::initialize() {
    log("Initializing GCMC simulation...");
    
    // Load parameters using setup module
    if (!setup_->loadParameters(config_.inputFile, *params_)) {
        log("Failed to load parameters");
        return false;
    }
    
    // Setup system
    if (!setup_->setupSystem(*params_, *state_)) {
        log("Failed to setup system");
        return false;
    }
    
    // Setup fragments
    if (!setup_->setupFragments(*params_, fragmentConfigs_, *reservoir_)) {
        log("Failed to setup fragments");
        return false;
    }
    
    // Setup acceptance calculator
    if (!setup_->setupAcceptance(*acceptance_, *params_)) {
        log("Failed to setup acceptance calculator");
        return false;
    }
    
    // Setup engine
    if (!setup_->setupEngine(*engine_, state_.get(), reservoir_.get(), *params_)) {
        log("Failed to setup engine");
        return false;
    }
    
    // Configure engine
    engine_->setAcceptanceCalculator(acceptance_.get());
    
    // Initialize core module
    core_->initialize(state_.get(), reservoir_.get(), 
                     engine_.get(), acceptance_.get());
    
    // Set random seed
    if (config_.randomSeed < 0) {
        std::random_device rd;
        config_.randomSeed = rd();
    }
    core_->setSeed(config_.randomSeed);
    
    // Configure core module
    core_->setVerbose(config_.verbose);
    
    // Validate setup
    if (!setup_->validateSetup(*params_, *state_)) {
        log("Setup validation failed");
        return false;
    }
    
    initialized_ = true;
    log("Initialization complete");
    
    // Print initial statistics
    if (config_.verbose) {
        printStatistics();
    }
    
    return true;
}

bool GCMCSimulationModular::run() {
    if (!initialized_) {
        log("Cannot run: simulation not initialized");
        return false;
    }
    
    log("Starting GCMC simulation...");
    running_ = true;
    
    // Get number of steps from parameters
    int numSteps = params_->get_mc_info().mc_steps;
    
    // Run simulation using core module
    bool success = core_->runSimulation(numSteps);
    
    running_ = false;
    
    if (success) {
        log("Simulation completed successfully");
    } else {
        log("Simulation failed");
    }
    
    return success;
}

void GCMCSimulationModular::finalize() {
    log("Finalizing simulation...");
    
    // Write final results
    // Write final results - simplified for now
    std::cout << "Simulation finalized" << std::endl;
    
    // Save final trajectory
    trajectoryWriter_->writeTrajectory(*state_, params_->get_mc_info().mc_steps, 
                        config_.outputPrefix + "_final.pdb");
    
    // Save final checkpoint
    checkpointManager_->saveCheckpoint(*state_, core_->getStatistics(), 
                       params_->get_mc_info().mc_steps,
                       config_.outputPrefix + "_final.chk");
    
    log("Simulation finalized");
}

void GCMCSimulationModular::stop() {
    if (running_) {
        running_ = false;
        log("Simulation stopped");
    }
}

bool GCMCSimulationModular::isRunning() const {
    return running_;
}

const StatisticsTracker& GCMCSimulationModular::getStatistics() const {
    return core_->getStatistics();
}

void GCMCSimulationModular::printStatistics() const {
    core_->getStatistics().printSummary(
        params_->get_mc_info().current_step);
}

void GCMCSimulationModular::writeTrajectory(int step) {
    trajectoryWriter_->writeTrajectory(*state_, step);
}

void GCMCSimulationModular::writeCheckpoint(int step) {
    checkpointManager_->saveCheckpoint(*state_, core_->getStatistics(), step);
}

bool GCMCSimulationModular::loadCheckpoint(const std::string& filename) {
    int step;
    StatisticsTracker stats;
    
    if (checkpointManager_->loadCheckpoint(*state_, stats, step, filename)) {
        params_->get_mc_info().current_step = step;
        // Note: Would need to transfer stats to core module
        log("Checkpoint loaded from step " + std::to_string(step));
        return true;
    }
    
    return false;
}

void GCMCSimulationModular::updateConfig(const Config& config) {
    config_ = config;
    
    // Update IO configuration
    io::output::TrajectoryWriter::Config trajConfig;
    trajConfig.prefix = config.outputPrefix;
    trajectoryWriter_->setConfig(trajConfig);
    
    io::output::CheckpointManager::Config chkConfig;
    chkConfig.prefix = config.outputPrefix;
    checkpointManager_->setConfig(chkConfig);
    
    // Update core configuration
    core_->setVerbose(config.verbose);
}

void GCMCSimulationModular::log(const std::string& message) const {
    if (config_.verbose) {
        std::cout << "[GCMCSimulation] " << message << std::endl;
    }
    
    // Log to file - simplified for now
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc