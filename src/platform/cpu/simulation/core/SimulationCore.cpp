#include "SimulationCore.hpp"
#include "../../energy/EnergyAPI.hpp"
#include "../../movement/MovementAPI.hpp"
#include "../../../../system/common/SystemLogger.hpp"
#include <stdexcept>
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace core {

using namespace pygcmc::system::common;
using namespace pygcmc::model::montecarlo;

SimulationCore::SimulationCore(const Config& config)
    : config_(config)
    , state_(nullptr)
    , engine_(nullptr)
    , initialized_(false)
    , running_(false)
    , stopRequested_(false)
    , currentStep_(0)
    , totalMoves_(0)
    , acceptedMoves_(0) {

    // Initialize random seed if specified
    if (config.randomSeed > 0) {
        // Set random seed through engine when initialized
    }
}

SimulationCore::~SimulationCore() {
    if (running_) {
        stop();
    }
}

void SimulationCore::initialize(MCState* state, movement::gcmc::GCMCEngine* engine) {
    if (!state || !engine) {
        throw std::invalid_argument("SimulationCore: state and engine must not be null");
    }

    state_ = state;
    engine_ = engine;
    initialized_ = true;

    // Apply energy method configuration
    applyEnergyMethod();

    // Reset statistics
    currentStep_ = 0;
    totalMoves_ = 0;
    acceptedMoves_ = 0;

    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("SimulationCore initialized with ",
                           config_.numSteps, " steps");
    }
}

bool SimulationCore::run(int nSteps) {
    checkInitialized();

    if (running_) {
        if (SystemLogger::isDebugEnabled()) {
            SystemLogger::debug("SimulationCore: Simulation already running");
        }
        return false;
    }

    running_ = true;
    stopRequested_ = false;

    int stepsToRun = (nSteps > 0) ? nSteps : config_.numSteps;
    int endStep = currentStep_ + stepsToRun;

    SystemLogger::info("Starting simulation for ", stepsToRun, " steps");

    // Equilibration phase
    if (currentStep_ < config_.equilibrationSteps) {
        int equilibrationSteps = std::min(config_.equilibrationSteps - currentStep_,
                                         stepsToRun);
        SystemLogger::info("Running equilibration for ", equilibrationSteps, " steps");
    }

    // Main simulation loop
    while (currentStep_ < endStep && !stopRequested_) {
        bool success = performStep();

        if (!success) {
            SystemLogger::error("SimulationCore: Step ", currentStep_, " failed");
            running_ = false;
            return false;
        }

        currentStep_++;

        // Call step callback if set
        if (stepCallback_ && currentStep_ % config_.statisticsInterval == 0) {
            stepCallback_(currentStep_, *state_);
        }
    }

    running_ = false;

    SystemLogger::info("Simulation completed. Total steps: ", currentStep_,
                      ", Acceptance rate: ", getAcceptanceRate());

    return true;
}

bool SimulationCore::performStep() {
    checkInitialized();

    // Perform a move
    bool moveAccepted = performMove();

    // Update statistics
    totalMoves_++;
    if (moveAccepted) {
        acceptedMoves_++;
    }

    return true;
}

bool SimulationCore::performMove() {
    checkInitialized();

    // Perform move through engine
    // TODO: GCMCEngine needs a performMove method that returns MovementResult
    // For now, just simulate a move
    platform::cpu::movement::MovementResult result;
    result.accepted = false;
    result.moveType = "translation";

    // Call move callback if set
    if (moveCallback_) {
        moveCallback_(result);
    }

    // Update internal statistics
    updateStatistics(result);

    return result.accepted;
}

void SimulationCore::stop() {
    stopRequested_ = true;
    running_ = false;

    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("SimulationCore: Stop requested");
    }
}

bool SimulationCore::isRunning() const {
    return running_;
}

void SimulationCore::reset() {
    currentStep_ = 0;
    totalMoves_ = 0;
    acceptedMoves_ = 0;
    stopRequested_ = false;
    running_ = false;
}

void SimulationCore::setStepCallback(StepCallback callback) {
    stepCallback_ = callback;
}

void SimulationCore::setMoveCallback(MoveCallback callback) {
    moveCallback_ = callback;
}

int SimulationCore::getCurrentStep() const {
    return currentStep_;
}

int SimulationCore::getTotalMoves() const {
    return totalMoves_;
}

int SimulationCore::getAcceptedMoves() const {
    return acceptedMoves_;
}

double SimulationCore::getAcceptanceRate() const {
    if (totalMoves_ == 0) {
        return 0.0;
    }
    return static_cast<double>(acceptedMoves_) / totalMoves_;
}

void SimulationCore::computeSystemEnergy() {
    checkInitialized();
    energy::computeSystemEnergy(*state_);
}

void SimulationCore::computeMovementEnergy() {
    checkInitialized();

    if (config_.energyMethod == EnergyMethod::DIRECT && config_.cutoff > 0) {
        energy::computeMovementEnergyCutoff(*state_);
    } else {
        energy::computeMovementEnergy(*state_);
    }
}

std::pair<double, double> SimulationCore::getEnergyComponents() const {
    checkInitialized();
    return energy::getTotalEnergyComponents(*state_);
}

void SimulationCore::updateConfig(const Config& config) {
    config_ = config;

    // Update random seed if changed
    if (config.randomSeed > 0) {
        // Set random seed through engine when available
    }

    // Apply new energy method
    if (initialized_) {
        applyEnergyMethod();
    }
}

void SimulationCore::checkInitialized() const {
    if (!initialized_) {
        throw std::runtime_error("SimulationCore: Not initialized");
    }
}

void SimulationCore::updateStatistics(const platform::cpu::movement::MovementResult& /* result */) {
    // Statistics are updated internally
    // Additional statistics tracking can be added here if needed
}

void SimulationCore::applyEnergyMethod() {
    // Set up energy calculation method based on configuration
    if (state_) {
        state_->info.cutoff = static_cast<float>(config_.cutoff / 10.0); // Convert to nm
    }

    // Additional energy method setup can be added here
    // For example, setting up Ewald, PME, or PGP parameters
}

} // namespace core
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc
