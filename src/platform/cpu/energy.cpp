// src/platform/cpu/energy.cpp
#include "energy.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

float computeNaiveNonbondedEnergy(model::MCState& state) {
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;

    // Reset energies for all residues
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        residue.energy_elec = 0.0f;
    }

    // Iterate through all active movement molecules
    for (const auto& movementInfo : state.movementResidues) {
        int start = movementInfo.startIndex;
        int count = movementInfo.activeCount;
        
        // For each active movement residue
        for (int i = start; i < start + count; ++i) {
            if (!residues[i].active) continue;
            
            // For each atom in the movement residue
            for (int atom_i = residues[i].atomStart; atom_i < residues[i].atomStart + residues[i].atomCount; ++atom_i) {
                int moveType = atoms[atom_i].type;
                
                // Compute interaction with atoms in all other active residues
                for (size_t j = 0; j < residues.size(); ++j) {
                    if (!residues[j].active || static_cast<int>(j) == i) continue;
                    
                    // For each atom in the other residue
                    for (int atom_j = residues[j].atomStart; atom_j < residues[j].atomStart + residues[j].atomCount; ++atom_j) {
                        int resType = atoms[atom_j].type;
                        
                        // Calculate distance between atoms
                        float dx = atoms[atom_j].x - atoms[atom_i].x;
                        float dy = atoms[atom_j].y - atoms[atom_i].y;
                        float dz = atoms[atom_j].z - atoms[atom_i].z;
                        float r = std::sqrt(dx*dx + dy*dy + dz*dz);
                        
                        size_t index = static_cast<size_t>(moveType * forcefield.maxTypes + resType);
                        if (index >= forcefield.ljEps.size()) continue;  // Safety check
                        
                        float eps = forcefield.ljEps[index];
                        float sigma = forcefield.ljSigma[index];
                        
                        // Debug output
                        printf("Computing energy between atoms %d (type %d) and %d (type %d)\n", 
                               atom_i, moveType, atom_j, resType);
                        printf("Distance r = %f\n", r);
                        printf("Using force field parameters: eps = %f, sigma = %f at index %zu\n", 
                               eps, sigma, index);
                        
                        // Calculate vdw energy: V = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
                        float term6 = std::pow(sigma / r, 6);
                        float term12 = term6 * term6;
                        float vdw_energy = eps * (term12 - 2.0f * term6);
                        
                        // Calculate electrostatic energy: V = q1*q2/(eps*r)
                        float q1 = atoms[atom_i].charge;
                        float q2 = atoms[atom_j].charge;
                        float elec_energy = q1 * q2 / r;  // Simple Coulomb
                        
                        printf("Computed vdw energy: %f\n", vdw_energy);
                        printf("Computed elec energy: %f\n", elec_energy);
                        
                        // Add the energies to the movement residue
                        residues[i].energy_vdw += vdw_energy;
                        residues[i].energy_elec += elec_energy;
                    }
                }
            }
        }
    }
    
    // Calculate total energy (divide by 2 because each interaction is counted twice)
    float totalEnergy = 0.0f;
    for (const auto& residue : residues) {
        totalEnergy += residue.energy_vdw + residue.energy_elec;
    }
    totalEnergy *= 0.5f;
    
    printf("Total energy: %f\n", totalEnergy);
    return totalEnergy;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc