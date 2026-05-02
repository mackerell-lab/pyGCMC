#pragma once

/**
 * @file DrudeMain.hpp
 * @brief Main header for Drude oscillator polarizable force field module
 *
 * This module implements the Drude oscillator model for molecular polarization,
 * following the same architectural patterns as other energy modules in PyGCMC.
 */

#include "DrudeInterface.hpp"
#include "DrudeCore.hpp"
#include "DrudeStructures.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Main entry point for Drude force calculations
 *
 * This class follows the same pattern as PMEComplete, PGPComplete, etc.
 * It provides a high-level interface for Drude oscillator calculations.
 */
class DrudeComplete {
public:
    /**
     * @brief Calculate Drude energy with SCF optimization
     * @param state Molecular state
     * @return Total Drude energy (harmonic + Thole corrections)
     */
    static double calculateEnergy(model::MCState& state);

    /**
     * @brief Calculate Drude energy using specific algorithm
     * @param state Molecular state
     * @param algorithm Optimization algorithm to use
     * @return Total Drude energy
     */
    static double calculateEnergy(model::MCState& state, DrudeAlgorithm algorithm);

    /**
     * @brief Set global Drude parameters
     * @param params SCF convergence parameters
     */
    static void setParameters(const DrudeSCFParams& params);

    /**
     * @brief Add a Drude particle to the system
     * @param particle Drude particle definition
     */
    static void addParticle(const DrudeParticle& particle);

    /**
     * @brief Add a Thole-screened pair interaction
     * @param pair Screened pair definition
     */
    static void addScreenedPair(const ScreenedPair& pair);

    /**
     * @brief Clear all Drude particles and pairs
     */
    static void clear();

    /**
     * @brief Get the number of Drude particles
     * @return Number of Drude particles
     */
    static size_t getNumParticles();

    /**
     * @brief Enable/disable ASPC history prediction
     * @param enable True to enable ASPC
     */
    static void enableASPC(bool enable);

    /**
     * @brief Check if ASPC is enabled
     * @return True if ASPC is enabled
     */
    static bool isASPCEnabled();

    /**
     * @brief Clear ASPC history
     */
    static void clearHistory();

    /**
     * @brief Get FastFBP optimizer instance for configuration
     * @return Pointer to FastFBP optimizer (nullptr if not available)
     */
    static DrudeFastFBP* getFastFBPOptimizer();

    /**
     * @brief Get Hybrid optimizer instance for configuration
     * @return Pointer to Hybrid optimizer (nullptr if not available)
     */
    static DrudeHybrid* getHybridOptimizer();

    /**
     * @brief Get MultiStage optimizer instance for configuration
     * @return Pointer to MultiStage optimizer (nullptr if not available)
     */
    static DrudeMultiStage* getMultiStageOptimizer();

    /**
     * @brief Get the DrudeCore instance
     * @return Reference to DrudeCore
     */
    static DrudeCore& getDrudeCore();
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
