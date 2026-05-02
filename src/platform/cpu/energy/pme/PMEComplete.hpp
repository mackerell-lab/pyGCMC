#pragma once

#include "model/ModelModule.hpp"
#include <memory>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace energy {
namespace pme {

/**
 * @brief Energy components structure
 */
struct EnergyComponents {
    double real{0.0};
    double reciprocal{0.0};
    double self{0.0};
    double vdw{0.0};
    double total{0.0};
};

/**
 * @brief Complete PME energy calculation with proper integration of all components
 *
 * This class provides a complete energy calculation
 * that includes both intermolecular and intramolecular interactions,
 * matching the behavior of reference implementations.
 */
class PMEComplete {
public:
    PMEComplete() = default;
    ~PMEComplete() = default;

    /**
     * @brief Compute complete system energy with all interactions
     *
     * This method calculates:
     * - Electrostatic energy (real + reciprocal + self)
     * - Van der Waals energy including intramolecular interactions
     * - Proper handling of exclusions and 1-4 interactions
     *
     * @param state The molecular state
     * @return Total energy struct with all components
     */
    EnergyComponents computeCompleteEnergy(model::MCState& state);

    /**
     * @brief Compute complete energy with cutoff method
     *
     * Similar to computeCompleteEnergy but uses cutoff instead of PME
     * for electrostatics. Useful for comparison and validation.
     *
     * @param state The molecular state
     * @return Total energy struct with all components
     */
    EnergyComponents computeCompleteEnergyCutoff(model::MCState& state);

private:
    /**
     * @brief Calculate intramolecular LJ interactions
     *
     * Computes LJ interactions within residues that are normally excluded
     * in the Fixed methods but needed for complete energy matching.
     *
     * @param state The molecular state
     * @return Intramolecular LJ energy
     */
    double calculateIntramolecularLJ(model::MCState& state);

    /**
     * @brief Apply proper exclusion rules for electrostatics
     *
     * Handles 1-2, 1-3, and 1-4 exclusions based on force field rules
     *
     * @param state The molecular state
     * @param energy Current energy components
     */
    void applyElectrostaticExclusions(model::MCState& state, EnergyComponents& energy);
};

} // namespace pme
} // namespace energy
} // namespace cpu
} // namespace platform
} // namespace pygcmc
