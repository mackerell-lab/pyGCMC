// src/platform/cpu/energy.cpp
#include "energy.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace pygcmc {
namespace platform {
namespace cpu {

// Debug flag to control output
static bool debug_output = false;

void computeNaiveNonbondedEnergy(model::MCState& state) {
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;

    // Validate force field setup
    size_t expected_size = static_cast<size_t>(forcefield.numMovementTypes) * 
                          static_cast<size_t>(forcefield.maxTypes);
    if (forcefield.ljEps.size() != expected_size) {
        std::stringstream ss;
        ss << "Force field parameters array size mismatch. Expected size "
           << expected_size
           << " (numMovementTypes * maxTypes), but got " << forcefield.ljEps.size();
        throw std::runtime_error(ss.str());
    }

    // Reset energies for all residues
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        residue.energy_elec = 0.0f;
    }

    // Iterate through all active movement molecules
    for (const auto& movementInfo : state.movementResidues) {
        size_t start = static_cast<size_t>(movementInfo.startIndex);
        size_t count = static_cast<size_t>(movementInfo.activeCount);
        
        // For each active movement residue
        for (size_t i = start; i < start + count; ++i) {
            if (!residues[i].active) continue;
            
            // For each atom in the movement residue
            for (size_t atom_i = static_cast<size_t>(residues[i].atomStart); 
                 atom_i < static_cast<size_t>(residues[i].atomStart + residues[i].atomCount); 
                 ++atom_i) {
                size_t moveType = static_cast<size_t>(atoms[atom_i].type);
                
                // Compute interaction with atoms in all other active residues
                for (size_t j = 0; j < residues.size(); ++j) {
                    if (!residues[j].active || j == i) continue;
                    
                    // For each atom in the other residue
                    for (size_t atom_j = static_cast<size_t>(residues[j].atomStart);
                         atom_j < static_cast<size_t>(residues[j].atomStart + residues[j].atomCount);
                         ++atom_j) {
                        size_t resType = static_cast<size_t>(atoms[atom_j].type);
                        
                        // Calculate distance between atoms
                        float dx = atoms[atom_j].x - atoms[atom_i].x;
                        float dy = atoms[atom_j].y - atoms[atom_i].y;
                        float dz = atoms[atom_j].z - atoms[atom_i].z;
                        float r = std::sqrt(dx*dx + dy*dy + dz*dz);
                        
                        size_t index = moveType * static_cast<size_t>(forcefield.maxTypes) + resType;
                        if (index >= forcefield.ljEps.size()) {
                            std::stringstream ss;
                            ss << "Invalid force field parameter index " << index 
                               << " for atom types " << moveType << " and " << resType
                               << ". This indicates a mismatch between atom types and force field parameters.";
                            throw std::runtime_error(ss.str());
                        }
                        
                        float eps = forcefield.ljEps[index];
                        float sigma = forcefield.ljSigma[index];
                        
                        platform::log(LogLevel::DEBUG, 
                            "Computing energy between atoms ", atom_i, " (type ", moveType, 
                            ") and ", atom_j, " (type ", resType, ")\n",
                            "Distance r = ", r, "\n",
                            "Using force field parameters: eps = ", eps, 
                            ", sigma = ", sigma, " at index ", index);
                        
                        // Calculate vdw energy: V = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
                        float term6 = std::pow(sigma / r, 6);
                        float term12 = term6 * term6;
                        float vdw_energy = eps * (term12 - 2.0f * term6);
                        
                        // Calculate electrostatic energy: V = q1*q2/(eps*r)
                        float q1 = atoms[atom_i].charge;
                        float q2 = atoms[atom_j].charge;
                        float elec_energy = q1 * q2 / r;  // Simple Coulomb
                        
                        platform::log(LogLevel::DEBUG, 
                            "Computed vdw energy: ", vdw_energy, "\n",
                            "Computed elec energy: ", elec_energy);
                        
                        // Add energies only to the movement residue
                        residues[i].energy_vdw += vdw_energy;
                        residues[i].energy_elec += elec_energy;
                    }
                }
            }
        }
    }
}

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable) {
    debug_output = enable;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc

