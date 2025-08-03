/**
 * @file DrudeMultiStage.cpp
 * @brief Implementation of multi-stage optimization: Direct → FastFBP → TCG → SCF
 */

#include "DrudeMultiStage.hpp"
#include "../common/EnergyConstants.hpp"
#include <chrono>
#include <cmath>
#include <numeric>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

DrudeMultiStage::DrudeMultiStage() {
    direct_ = std::make_unique<DrudeDirectPolarization>();
    fastFBP_ = std::make_unique<DrudeFastFBP>();
    tcg_ = std::make_unique<DrudeTCG>();
    scf_ = std::make_unique<DrudeSCF>();
    
    // Configure TCG based on literature recommendations
    tcg_->setIterations(3);  // TCG-3 optimal for speed/accuracy
    tcg_->setUseChebyshev(true);
}

bool DrudeMultiStage::optimize(
    model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const DrudeSCFParams& params
) {
    // auto start_total = std::chrono::high_resolution_clock::now();
    
    // Initialize statistics
    stats_ = MultiStageStats{};
    
    // Analyze system characteristics
    analyzeSystem(state, particles);
    
    // Adapt configuration based on system
    if (config_.adaptiveMode) {
        adaptConfiguration(state, particles);
    }
    
    // Stage 1: Direct Polarization
    auto start_direct = std::chrono::high_resolution_clock::now();
    
    if (useDirectStage_) {
        direct_->optimize(state, particles, screenedPairs, params);
        
        auto end_direct = std::chrono::high_resolution_clock::now();
        stats_.directTime = std::chrono::duration<double>(end_direct - start_direct).count();
        
        // Calculate error after direct stage
        stats_.errorAfterDirect = calculateError(state, particles);
        
        // If system is very sparse and direct is already good, might skip other stages
        if (stats_.systemDensity < 0.1 && stats_.errorAfterDirect < 0.001) {
            stats_.converged = true;
            return true;
        }
    }
    
    // Stage 2: FastFBP
    auto start_fbp = std::chrono::high_resolution_clock::now();
    
    if (useFBPStage_) {
        // Configure FastFBP
        fastFBP_->setIterationMode(DrudeFastFBP::IterationMode::Fixed);
        
        // Determine optimal number of iterations
        int fbpIterations = determineFBPIterations(
            stats_.systemDensity,
            stats_.avgPolarizability,
            particles.size()
        );
        
        // For dynamic convergence checking
        double previousError = stats_.errorAfterDirect;
        
        // Run FastFBP iterations with monitoring
        for (int iter = 0; iter < fbpIterations; ++iter) {
            fastFBP_->setIterations(1);
            fastFBP_->optimize(state, particles, screenedPairs, params);
            
            double currentError = calculateError(state, particles);
            double convergenceRate = (previousError - currentError) / previousError;
            
            stats_.fbpIterations = iter + 1;
            
            // Check if we should switch to SCF early
            if (shouldSwitchToSCF(currentError, convergenceRate, iter + 1)) {
                break;
            }
            
            previousError = currentError;
        }
        
        auto end_fbp = std::chrono::high_resolution_clock::now();
        stats_.fbpTime = std::chrono::duration<double>(end_fbp - start_fbp).count();
        stats_.errorAfterFBP = previousError;
        
        // Check if we should skip TCG and go directly to SCF
        if (!config_.enableTCG || stats_.errorAfterFBP > config_.tcgErrorThreshold) {
            // Skip TCG if error is still too large or TCG disabled
            useTCGStage_ = false;
        }
    }
    
    // Stage 3: TCG Refinement (if enabled and appropriate)
    auto start_tcg = std::chrono::high_resolution_clock::now();
    
    if (useTCGStage_ && config_.enableTCG) {
        // Configure TCG iterations
        tcg_->setIterations(config_.tcgIterations);
        
        // Run TCG
        DrudeSCFParams tcgParams = params;
        tcgParams.maxIterations = config_.tcgIterations;  // Fixed iterations
        
        tcg_->optimize(state, particles, screenedPairs, tcgParams);
        
        auto end_tcg = std::chrono::high_resolution_clock::now();
        stats_.tcgTime = std::chrono::duration<double>(end_tcg - start_tcg).count();
        stats_.tcgIterations = config_.tcgIterations;
        
        // Calculate error after TCG
        stats_.errorAfterTCG = calculateError(state, particles);
        
        // Check if TCG is sufficient
        if (stats_.errorAfterTCG < params.tolerance && !config_.requireSCF) {
            stats_.finalError = stats_.errorAfterTCG;
            stats_.converged = true;
            return true;
        }
    } else {
        stats_.errorAfterTCG = stats_.errorAfterFBP;
        stats_.tcgTime = 0.0;
        stats_.tcgIterations = 0;
    }
    
    // Stage 4: SCF Fine-tuning (if required)
    auto start_scf = std::chrono::high_resolution_clock::now();
    
    if (!useSCFStage_ || (!config_.requireSCF && stats_.errorAfterTCG < params.tolerance)) {
        // Skip SCF if not required
        stats_.finalError = stats_.errorAfterTCG;
        stats_.converged = true;
        stats_.scfTime = 0.0;
        stats_.scfIterations = 0;
        return true;
    }
    
    // Create modified SCF parameters
    DrudeSCFParams scfParams = params;
    scfParams.maxIterations = config_.maxSCFIterations;
    
    // If we're already close, can use looser tolerance
    double currentError = useTCGStage_ ? stats_.errorAfterTCG : stats_.errorAfterFBP;
    if (currentError < 0.01) {
        scfParams.tolerance = std::max(config_.scfTolerance, currentError * 0.1);
    }
    
    // Run SCF
    bool converged = scf_->optimize(state, particles, screenedPairs, scfParams);
    
    auto end_scf = std::chrono::high_resolution_clock::now();
    stats_.scfTime = std::chrono::duration<double>(end_scf - start_scf).count();
    stats_.scfIterations = scf_->getIterationCount();
    
    // Final error
    stats_.finalError = calculateError(state, particles);
    stats_.converged = converged;
    
    // Debug output (if enabled)
    if (false) {  // Set to true for debugging
        std::cout << "Multi-stage optimization complete:" << std::endl;
        std::cout << "  System density: " << stats_.systemDensity << " g/cm³" << std::endl;
        std::cout << "  Avg polarizability: " << stats_.avgPolarizability << " nm³" << std::endl;
        std::cout << "  Direct: " << stats_.directTime*1000 << " ms, error: " << stats_.errorAfterDirect << std::endl;
        std::cout << "  FastFBP: " << stats_.fbpTime*1000 << " ms (" << stats_.fbpIterations << " iter), error: " << stats_.errorAfterFBP << std::endl;
        if (stats_.tcgIterations > 0) {
            std::cout << "  TCG: " << stats_.tcgTime*1000 << " ms (" << stats_.tcgIterations << " iter), error: " << stats_.errorAfterTCG << std::endl;
        }
        std::cout << "  SCF: " << stats_.scfTime*1000 << " ms (" << stats_.scfIterations << " iter), error: " << stats_.finalError << std::endl;
        std::cout << "  Total time: " << (stats_.directTime + stats_.fbpTime + stats_.tcgTime + stats_.scfTime)*1000 << " ms" << std::endl;
    }
    
    return converged;
}

void DrudeMultiStage::analyzeSystem(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles
) {
    // Calculate system density
    double volume = state.info.box[0] * state.info.box[1] * state.info.box[2];  // nm³
    double mass = 0.0;
    
    // Rough estimate: assume average atomic mass of 15 amu
    mass = state.activeAtomCount * 15.0;  // amu
    
    // Convert to g/cm³: 1 amu/nm³ = 1.66054 g/cm³
    stats_.systemDensity = (mass / volume) * 1.66054;
    
    // Calculate average polarizability
    if (!particles.empty()) {
        double totalPolarizability = 0.0;
        for (const auto& p : particles) {
            totalPolarizability += p.polarizability;
        }
        stats_.avgPolarizability = totalPolarizability / particles.size();
    }
}

void DrudeMultiStage::adaptConfiguration(
    const model::MCState& /*state*/,
    const std::vector<DrudeParticle>& particles
) {
    // Adapt based on density
    if (stats_.systemDensity < 0.5) {
        // Low density: Direct is already good, minimal FastFBP
        config_.minFBPIterations = 2;
        config_.maxFBPIterations = 5;
        config_.fbpCutoffFactor = 0.6;
        config_.tcgIterations = 2;  // Fewer TCG iterations needed
        config_.requireSCF = false;  // Often don't need SCF
    } else if (stats_.systemDensity < 1.0) {
        // Medium density: Standard configuration
        config_.minFBPIterations = 3;
        config_.maxFBPIterations = 7;
        config_.fbpCutoffFactor = 0.8;
        config_.tcgIterations = 3;  // TCG-3 optimal
        config_.requireSCF = false;  // TCG usually sufficient
    } else {
        // High density: Need more iterations
        config_.minFBPIterations = 5;
        config_.maxFBPIterations = 10;
        config_.fbpCutoffFactor = 1.0;
        config_.tcgIterations = 4;  // More TCG iterations
        config_.requireSCF = true;   // May need SCF refinement
    }
    
    // Adapt based on polarizability
    if (stats_.avgPolarizability > 0.002) {
        // High polarizability: Need more iterations
        config_.maxFBPIterations += 2;
        config_.tcgIterations += 1;  // Extra TCG iteration
        config_.switchToSCFError *= 0.5;  // Switch earlier
        config_.requireSCF = true;  // Likely need SCF
    }
    
    // Adapt based on system size
    if (particles.size() > 1000) {
        // Large system: Relax tolerances slightly
        config_.scfTolerance *= 2;
        config_.switchToSCFError *= 2;
        config_.tcgErrorThreshold *= 2;  // Relax TCG threshold
        config_.tcgIterations = std::min(config_.tcgIterations, 3);  // Limit TCG
    }
}

double DrudeMultiStage::calculateError(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles
) const {
    // Calculate RMS displacement as error metric
    double totalError = 0.0;
    
    for (const auto& particle : particles) {
        const auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        
        // Apply PBC
        const auto& box = state.info.box;
        if (dx > box[0]/2) dx -= box[0];
        if (dx < -box[0]/2) dx += box[0];
        if (dy > box[1]/2) dy -= box[1];
        if (dy < -box[1]/2) dy += box[1];
        if (dz > box[2]/2) dz -= box[2];
        if (dz < -box[2]/2) dz += box[2];
        
        double displacement2 = dx*dx + dy*dy + dz*dz;
        totalError += displacement2;
    }
    
    return std::sqrt(totalError / particles.size());
}

int DrudeMultiStage::determineFBPIterations(
    double density,
    double polarizability,
    size_t systemSize
) const {
    // Base iterations from configuration
    int iterations = config_.minFBPIterations;
    
    // Adjust based on density
    if (density > 1.5) {
        iterations += 2;
    } else if (density > 1.0) {
        iterations += 1;
    }
    
    // Adjust based on polarizability
    if (polarizability > 0.002) {
        iterations += 2;
    } else if (polarizability > 0.001) {
        iterations += 1;
    }
    
    // Adjust based on system size
    if (systemSize > 500) {
        iterations += 1;
    }
    
    // Cap at maximum
    return std::min(iterations, config_.maxFBPIterations);
}

bool DrudeMultiStage::shouldSwitchToSCF(
    double currentError,
    double convergenceRate,
    int fbpIterations
) const {
    // Switch if error is already small enough for TCG
    if (currentError < config_.tcgErrorThreshold) {
        return true;
    }
    
    // Switch if convergence has stalled
    if (fbpIterations >= 3 && convergenceRate < config_.switchToSCFRate) {
        return true;
    }
    
    // Switch if reached minimum iterations
    if (fbpIterations >= config_.minFBPIterations && currentError < 2 * config_.switchToSCFError) {
        return true;
    }
    
    // Don't switch if haven't reached minimum iterations
    return false;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc