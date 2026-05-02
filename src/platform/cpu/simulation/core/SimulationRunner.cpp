#include "SimulationRunner.hpp"
#include "../impl/GCMCSimulation.hpp"
#include <iostream>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

SimulationRunner::SimulationRunner(GCMCSimulation* parent)
    : parent_(parent)
    , state_(nullptr)
    , reservoir_(nullptr)
    , engine_(nullptr)
    , acceptance_(nullptr)
    , uniform_(0.0, 1.0)
    , verbose_(false)
    , printFrequency_(1000)
    , trajectoryFrequency_(10000)
    , checkpointFrequency_(100000)
    , enableStatistics_(true)
    , statisticsInterval_(1000) {
}

SimulationRunner::~SimulationRunner() {
    // Components are not owned, so no deletion needed
}

void SimulationRunner::initialize(model::montecarlo::MCState* state,
                               movement::MultiTypeReservoir* reservoir,
                               movement::gcmc::GCMCEngine* engine,
                               movement::gcmc::GCMCAcceptance* acceptance) {
    state_ = state;
    reservoir_ = reservoir;
    engine_ = engine;
    acceptance_ = acceptance;
}

void SimulationRunner::setSeed(unsigned int seed) {
    rng_.seed(seed);
}

bool SimulationRunner::runSimulation(int numSteps) {
    if (!state_ || !reservoir_ || !engine_) {
        std::cerr << "SimulationCore: Not initialized" << std::endl;
        return false;
    }

    auto startTime = std::chrono::steady_clock::now();

    for (int step = 1; step <= numSteps; ++step) {
        // Perform MC step
        if (!performMCStep()) {
            std::cerr << "SimulationCore: Failed at step " << step << std::endl;
            return false;
        }

        // Print statistics
        if (step % printFrequency_ == 0 && step > 0) {
            statistics_.printSummary(step);
        }

        // Save trajectory
        if (step % trajectoryFrequency_ == 0) {
            // Skip for now - parent should handle this
        }

        // Save checkpoint
        if (step % checkpointFrequency_ == 0) {
            // Skip for now - parent should handle this
        }
    }

    auto endTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;

    if (verbose_) {
        std::cout << "Simulation completed in " << elapsed.count() << " seconds" << std::endl;
        std::cout << "Performance: " << numSteps / elapsed.count() << " steps/second" << std::endl;
    }

    return true;
}

bool SimulationRunner::performMCStep() {
    // Select move type
    MoveType moveType = selectMoveType();

    switch (moveType) {
        case INSERT:
            performInsertion(selectFragmentType());
            break;
        case DELETE:
            performDeletion();
            break;
        case TRANSLATE:
            performTranslation();
            break;
        case ROTATE:
            performRotation();
            break;
    }

    // Update energy statistics if enabled
    if (enableStatistics_ && statistics_.getTotalSteps() % statisticsInterval_ == 0) {
        updateEnergyStatistics();
    }

    return true;
}

bool SimulationRunner::performInsertion(int fragType) {
    if (fragType < 0) return false;

    auto result = engine_->attemptInsertion(fragType);

    // Get fragment name from parent
    std::string fragName = "fragment";
    if (parent_) {
        auto fragInfo = parent_->getFragmentInfo();
        for (const auto& frag : fragInfo) {
            if (frag.typeId == fragType) {
                fragName = frag.name;
                break;
            }
        }
    }

    statistics_.recordMove("insert", fragName, result.accepted);
    return result.accepted;
}

bool SimulationRunner::performDeletion() {
    if (reservoir_->getActiveCount() == 0) {
        return false;
    }

    int fragType = selectActiveFragment();
    if (fragType < 0) return false;

    auto result = engine_->attemptDeletion(fragType);

    // Get fragment name
    std::string fragName = "fragment";
    if (parent_) {
        auto fragInfo = parent_->getFragmentInfo();
        for (const auto& frag : fragInfo) {
            if (frag.typeId == fragType) {
                fragName = frag.name;
                break;
            }
        }
    }

    statistics_.recordMove("delete", fragName, result.accepted);
    return result.accepted;
}

bool SimulationRunner::performTranslation() {
    if (reservoir_->getActiveCount() == 0) {
        return false;
    }

    auto activeIndices = reservoir_->getActiveInstances();
    if (activeIndices.empty()) return false;

    int idx = activeIndices[rng_() % activeIndices.size()];
    auto result = engine_->attemptTranslation(idx);

    statistics_.recordMove("translate", "fragment", result.accepted);
    return result.accepted;
}

bool SimulationRunner::performRotation() {
    if (reservoir_->getActiveCount() == 0) {
        return false;
    }

    auto activeIndices = reservoir_->getActiveInstances();
    if (activeIndices.empty()) return false;

    int idx = activeIndices[rng_() % activeIndices.size()];
    auto result = engine_->attemptRotation(idx);

    statistics_.recordMove("rotate", "fragment", result.accepted);
    return result.accepted;
}

SimulationRunner::MoveType SimulationRunner::selectMoveType() {
    double r = uniform_(rng_);

    // Equal probability for now
    if (r < 0.25) return INSERT;
    else if (r < 0.50) return DELETE;
    else if (r < 0.75) return TRANSLATE;
    else return ROTATE;
}

int SimulationRunner::selectFragmentType() {
    if (!parent_) return 0;

    auto fragInfo = parent_->getFragmentInfo();
    if (fragInfo.empty()) return -1;

    // For now, select randomly with equal probability
    return fragInfo[rng_() % fragInfo.size()].typeId;
}

int SimulationRunner::selectActiveFragment() {
    if (!reservoir_) return -1;

    auto activeIndices = reservoir_->getActiveInstances();
    if (activeIndices.empty()) return -1;

    // Randomly select an active fragment
    // Note: This is simplified, in real implementation we'd track this properly
    return 0;  // Return type 0 for now
}

double SimulationRunner::calculateSystemEnergy() {
    if (!engine_) return 0.0;
    return engine_->calculateSystemEnergy();
}

void SimulationRunner::updateEnergyStatistics() {
    double energy = calculateSystemEnergy();
    statistics_.recordEnergy(energy);
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc
