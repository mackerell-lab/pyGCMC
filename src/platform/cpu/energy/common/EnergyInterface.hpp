#pragma once

#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Energy calculation method enumeration
enum class EnergyMethod {
    DIRECT,  // Direct calculation method
    EWALD,   // Ewald summation method
    PME      // Particle Mesh Ewald method
};

// Unified energy calculation interface
void computeSystemEnergy(model::MCState& state, 
                         EnergyMethod method = EnergyMethod::DIRECT,
                         bool use_cutoff = false, 
                         bool use_pbc = false);

void computeMovementEnergy(model::MCState& state, 
                          EnergyMethod method = EnergyMethod::DIRECT,
                          bool use_cutoff = false, 
                          bool use_pbc = false);

// Forward declarations for all direct calculation functions to avoid ambiguity
void computeMovementEnergy(model::MCState& state);
void computeMovementEnergyCutoff(model::MCState& state);
void computeSystemEnergy(model::MCState& state);
void computeSystemEnergyCutoff(model::MCState& state);
void computeSystemEnergyPBC(model::MCState& state);
void computeSystemEnergyPBCCutoff(model::MCState& state);
void computeSystemVdwEnergyCutoff(model::MCState& state);

// Declare Direct, Ewald and PME calculation functions (internal use)
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

// Note: Ewald, PME, PGP and debug functions are defined in their respective modules
// This interface only provides the direct calculation functions and unified interface

// Clear PME state function
void clearPMEState();

// Get total energy function
inline double getTotalEnergy(const model::MCState& state) {
    double total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            total += residue.energy_vdw + residue.energy_elec;
        }
    }
    return total;
}

inline double getEwaldTotalEnergy(const model::MCState& state) {
    // Calculate energy from all residues (including vdw energy and real-space electrostatic energy)
    double residue_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            residue_total += residue.energy_vdw + residue.energy_elec;
        }
    }
    
    // Add reciprocal-space energy and self-energy from EwaldEnergy
    return residue_total + state.ewald_energy.reciprocal + state.ewald_energy.self;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 