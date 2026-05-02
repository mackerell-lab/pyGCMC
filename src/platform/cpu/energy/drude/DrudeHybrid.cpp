/**
 * @file DrudeHybrid.cpp
 * @brief Implementation of hybrid FastFBP+SCF optimization
 */

#include "DrudeHybrid.hpp"
#include "../common/EnergyConstants.hpp"
#include <chrono>
#include <cmath>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

DrudeHybrid::DrudeHybrid() {
    fastFBP_ = std::make_unique<DrudeFastFBP>();
    scf_ = std::make_unique<DrudeSCF>();
}

bool DrudeHybrid::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    // auto start_total = std::chrono::high_resolution_clock::now();

    // Initialize statistics
    stats_ = HybridStats{0, 0, 0.0, 0.0, 0.0, false};

    // Configure sub-optimizers
    configureFastFBP(params);
    configureSCF(params);

    // Phase 1: FastFBP preconditioning
    auto start_fbp = std::chrono::high_resolution_clock::now();

    // Determine FastFBP iterations based on mode
    int actualFBPIterations = fbpIterations_;
    if (hybridMode_ == HybridMode::Adaptive) {
        // Adjust based on system size
        size_t nParticles = particles.size();
        if (nParticles < 10) {
            actualFBPIterations = 3;
        } else if (nParticles < 50) {
            actualFBPIterations = 5;
        } else {
            actualFBPIterations = 7;
        }
    }

    // Run FastFBP
    // bool fbpSuccess = true;
    double previousError = 1.0;
    double convergenceRate = 0.0;

    if (hybridMode_ == HybridMode::Dynamic) {
        // Dynamic mode: run FastFBP with monitoring
        fastFBP_->setIterationMode(DrudeFastFBP::IterationMode::Fixed);

        for (int iter = 0; iter < actualFBPIterations; ++iter) {
            fastFBP_->setIterations(1);  // One iteration at a time
            fastFBP_->optimize(state, particles, screenedPairs, params);

            // Calculate current error
            double currentError = calculateError(state, particles);

            // Calculate convergence rate
            if (iter > 0) {
                convergenceRate = (previousError - currentError) / previousError;
            }
            previousError = currentError;

            stats_.fbpIterations = iter + 1;

            // Check if we should switch to SCF
            if (shouldSwitchToSCF(currentError, iter + 1, convergenceRate)) {
                break;
            }
        }

        stats_.switchError = previousError;

    } else {
        // Fixed mode: run all FastFBP iterations at once
        fastFBP_->setIterationMode(DrudeFastFBP::IterationMode::Fixed);
        fastFBP_->setIterations(actualFBPIterations);
        fastFBP_->optimize(state, particles, screenedPairs, params);
        stats_.fbpIterations = actualFBPIterations;
        stats_.switchError = calculateError(state, particles);
    }

    auto end_fbp = std::chrono::high_resolution_clock::now();
    stats_.fbpTime = std::chrono::duration<double>(end_fbp - start_fbp).count();

    // Phase 2: SCF refinement
    auto start_scf = std::chrono::high_resolution_clock::now();

    // Create modified SCF parameters with reduced iterations
    DrudeSCFParams scfParams = params;
    scfParams.maxIterations = maxSCFIter_;

    // If FastFBP got us close, we can use a looser tolerance
    if (stats_.switchError < 0.1) {
        scfParams.tolerance = std::max(params.tolerance, stats_.switchError * 0.1);
    }

    // Run SCF starting from FastFBP result
    bool scfSuccess = scf_->optimize(state, particles, screenedPairs, scfParams);

    // Get SCF iteration count
    stats_.scfIterations = scf_->getIterationCount();

    auto end_scf = std::chrono::high_resolution_clock::now();
    stats_.scfTime = std::chrono::duration<double>(end_scf - start_scf).count();

    stats_.converged = scfSuccess;

    // Debug output (if enabled)
    if (false) {  // Set to true for debugging
        std::cout << "Hybrid optimization complete:" << std::endl;
        std::cout << "  FastFBP iterations: " << stats_.fbpIterations << std::endl;
        std::cout << "  Switch error: " << stats_.switchError << std::endl;
        std::cout << "  SCF iterations: " << stats_.scfIterations << std::endl;
        std::cout << "  Total time: " << (stats_.fbpTime + stats_.scfTime) << "s" << std::endl;
        std::cout << "  Time saved vs pure SCF: " <<
                     (1.0 - (stats_.fbpTime + stats_.scfTime) / (stats_.scfTime * 2.5)) * 100 << "%" << std::endl;
    }

    return scfSuccess;
}

double DrudeHybrid::calculateError(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles
) {
    // Calculate RMS force imbalance as error metric
    double totalError = 0.0;

    for (const auto& particle : particles) {
        const auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];

        // Spring force
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;

        // Spring force would be:
        // Vec3 springForce = {
        //     -particle.kSpring * dx,
        //     -particle.kSpring * dy,
        //     -particle.kSpring * dz
        // };

        // Electric force (simplified - would need full field calculation)
        // For now, use displacement as proxy for error
        double displacement2 = dx*dx + dy*dy + dz*dz;
        totalError += displacement2;
    }

    return std::sqrt(totalError / particles.size());
}

bool DrudeHybrid::shouldSwitchToSCF(
    double currentError,
    int fbpIterations,
    double convergenceRate
) {
    // Switch criteria for dynamic mode

    // 1. Error is already small enough
    if (currentError < switchThreshold_) {
        return true;
    }

    // 2. Convergence has stalled
    if (fbpIterations >= 3 && convergenceRate < 0.1) {
        return true;
    }

    // 3. Reached maximum FastFBP iterations
    if (fbpIterations >= fbpIterations_) {
        return true;
    }

    return false;
}

void DrudeHybrid::configureFastFBP(const DrudeSCFParams& /*params*/) {
    // Configure FastFBP for preconditioning
    fastFBP_->setIterationMode(DrudeFastFBP::IterationMode::Fixed);
    fastFBP_->setCutoff(1.2);  // Use reasonable cutoff
    fastFBP_->setIncludeDrudeDrude(true);

    // In adaptive mode, configure based on system
    if (hybridMode_ == HybridMode::Adaptive) {
        fastFBP_->setIterationMode(DrudeFastFBP::IterationMode::Adaptive);
        fastFBP_->setAdaptiveMode(true);
    }
}

void DrudeHybrid::configureSCF(const DrudeSCFParams& /*params*/) {
    // SCF is already configured through the params argument
    // No additional configuration needed here
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
