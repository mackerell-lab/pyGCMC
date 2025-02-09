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
                          static_cast<size_t>(forcefield.numTotalTypes);
    if (forcefield.ljEps.size() != expected_size || forcefield.ljSigma.size() != expected_size) {
        std::stringstream ss;
        ss << "Force field parameters array size mismatch. Expected size "
           << expected_size
           << " (numMovementTypes * numTotalTypes), but got eps=" << forcefield.ljEps.size()
           << " sigma=" << forcefield.ljSigma.size();
        throw std::runtime_error(ss.str());
    }

    // Reset energies for all residues
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        residue.energy_elec = 0.0f;
    }

    // Iterate through all movement molecule groups
    for (const auto& movementInfo : state.movementResidues) {
        // Calculate energies only for active movement residues
        for (int i = movementInfo.startIndex; 
             i < movementInfo.startIndex + movementInfo.activeCount; ++i) {
            if (!residues[i].active) continue;
            
            // For each atom in the movement residue
            for (int atom_i = residues[i].atomStart; 
                 atom_i < residues[i].atomStart + residues[i].atomCount; 
                 ++atom_i) {
                int moveType = atoms[atom_i].type;
                
                // Validate movement atom type
                if (moveType >= forcefield.numMovementTypes) {
                    std::stringstream ss;
                    ss << "Movement atom type " << moveType << " out of range. "
                       << "Maximum allowed type is " << (forcefield.numMovementTypes - 1)
                       << " for atom " << atom_i << " in residue " << i 
                       << " (" << movementInfo.resName << ")";
                    throw std::runtime_error(ss.str());
                }
                
                // Compute interaction with atoms in all other active residues
                for (int j = 0; j < state.activeResidueCount; ++j) {
                    // Skip if:
                    // 1. Same residue
                    // 2. Residue is inactive
                    if (j == i || !residues[j].active) continue;
                    
                    // For each atom in the other residue
                    for (int atom_j = residues[j].atomStart;
                         atom_j < residues[j].atomStart + residues[j].atomCount;
                         ++atom_j) {
                        int resType = atoms[atom_j].type;
                        
                        // Validate residue atom type
                        if (resType >= forcefield.numTotalTypes) {
                            std::stringstream ss;
                            ss << "Residue atom type " << resType << " out of range. "
                               << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
                               << " for atom " << atom_j << " in residue " << j;
                            throw std::runtime_error(ss.str());
                        }
                        
                        // Calculate distance between atoms
                        float dx = atoms[atom_j].x - atoms[atom_i].x;
                        float dy = atoms[atom_j].y - atoms[atom_i].y;
                        float dz = atoms[atom_j].z - atoms[atom_i].z;
                        float r2 = dx*dx + dy*dy + dz*dz;  // nm^2
                        float r = std::sqrt(r2);  // nm
                        
                        // Get force field parameters
                        int param_index = moveType * forcefield.numTotalTypes + resType;
                        float eps = forcefield.ljEps[param_index];
                        float sigma = forcefield.ljSigma[param_index];
                        
                        // Calculate vdw energy: V = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
                        float sigma_r = sigma / r;
                        float term6 = std::pow(sigma_r, 6);
                        float term12 = term6 * term6;
                        float vdw_energy = eps * (term12 - 2.0f * term6);
                        
                        // Calculate electrostatic energy: V = q1*q2/r
                        float q1 = atoms[atom_i].charge;
                        float q2 = atoms[atom_j].charge;
                        float elec_energy = q1 * q2 / r;
                        
                        // Add energies to the movement residue
                        residues[i].energy_vdw += vdw_energy;
                        residues[i].energy_elec += elec_energy;
                        
                        if (debug_output) {
                            platform::log(LogLevel::DEBUG,
                                "Interaction between residues ", i, "(", movementInfo.resName, ") and ", j, "\n",
                                "  Atoms: ", atom_i, "(type ", moveType, ") - ", 
                                atom_j, "(type ", resType, ")\n",
                                "  Distance: ", r, " nm\n",
                                "  Parameters: eps=", eps, " sigma=", sigma, "\n",
                                "  Energies: vdw=", vdw_energy, " elec=", elec_energy);
                        }
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

