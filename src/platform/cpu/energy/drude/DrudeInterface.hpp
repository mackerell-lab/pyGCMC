#pragma once

/**
 * @file DrudeInterface.hpp
 * @brief Interface definitions for Drude force calculations
 */

#include "model/ModelModule.hpp"
#include "../common/EnergyInterface.hpp"
#include <vector>
#include <array>

namespace pygcmc {
namespace platform {
namespace cpu {

// Forward declarations
struct DrudeParticle;
struct ScreenedPair;
struct DrudeSCFParams;
enum class DrudeAlgorithm;

// Type alias for 3D vectors
using Vec3 = std::array<double, 3>;

/**
 * @brief Abstract interface for Drude optimization algorithms
 * 
 * This follows the strategy pattern to allow different optimization
 * algorithms (SCF, OPT3, FBP) to be used interchangeably.
 */
class DrudeOptimizer {
public:
    
    virtual ~DrudeOptimizer() = default;
    
    /**
     * @brief Optimize Drude particle positions
     * @param state Molecular state (will be modified)
     * @param particles Vector of Drude particles
     * @param screenedPairs Vector of Thole-screened pairs
     * @param params SCF parameters
     * @return true if converged, false otherwise
     */
    virtual bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) = 0;
    
    /**
     * @brief Get algorithm name for debugging
     * @return Algorithm name
     */
    virtual const char* getName() const = 0;
};

/**
 * @brief Interface for Drude force calculations
 * 
 * This class manages Drude particles and delegates optimization
 * to specific algorithm implementations.
 */
class DrudeInterface {
public:
    virtual ~DrudeInterface() = default;
    
    /**
     * @brief Calculate Drude energy after optimization
     * @param state Molecular state
     * @return Total Drude energy (harmonic + Thole)
     */
    virtual double calculateEnergy(model::MCState& state) = 0;
    
    /**
     * @brief Calculate forces from Drude interactions
     * @param state Molecular state
     * @param forces Force vector to update
     */
    virtual void calculateForces(
        model::MCState& state,
        std::vector<Vec3>& forces
    ) = 0;
    
    /**
     * @brief Add a Drude particle
     * @param particle Drude particle parameters
     * @return Index of added particle
     */
    virtual int addParticle(const DrudeParticle& particle) = 0;
    
    /**
     * @brief Add a Thole-screened pair
     * @param pair Screened pair parameters
     */
    virtual void addScreenedPair(const ScreenedPair& pair) = 0;
    
    /**
     * @brief Set optimization algorithm
     * @param algorithm Algorithm to use
     */
    virtual void setAlgorithm(DrudeAlgorithm algorithm) = 0;
    
    /**
     * @brief Set SCF parameters
     * @param params Convergence parameters
     */
    virtual void setParameters(const DrudeSCFParams& params) = 0;
    
    /**
     * @brief Clear all particles and pairs
     */
    virtual void clear() = 0;
    
    /**
     * @brief Get number of Drude particles
     * @return Number of particles
     */
    virtual size_t getNumParticles() const = 0;
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc