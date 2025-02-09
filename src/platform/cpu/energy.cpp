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

/**
 * @brief Coulomb constant in GROMACS MD units [kJ·nm/mol/e²]
 * 
 * k_c = 1/(4*π*ε₀) = 138.935458 kJ·nm/mol/e²
 * 
 * Unit analysis:
 * - ε₀ (vacuum permittivity) = 8.8541878128e-12 C²/(J·m)
 * - 1 kJ = 1000 J
 * - 1 nm = 1e-9 m
 * - 1 e = 1.60217663e-19 C
 * - N_A (Avogadro constant) = 6.02214076e23 mol⁻¹
 */
const float COULOMB = 138.935458f;

/**
 * @brief Safety parameters for energy calculation
 * 
 * !!! CRITICAL: Distance handling for Monte Carlo simulation !!!
 * 
 * MIN_SAFE_DISTANCE: Minimum allowed distance (1% of sigma)
 * - !!! Prevents numerical instability and infinity at r = 0
 * - !!! Essential for Monte Carlo sampling near contact
 * - !!! Implements soft core potential for r < MIN_SAFE_DISTANCE
 * 
 * MAX_SAFE_ENERGY: Maximum allowed energy per interaction
 * - !!! Prevents numerical overflow in Metropolis criterion
 * - !!! Keeps energies finite for stable MC sampling
 * - !!! Especially important for Coulomb interactions at small r
 */
const float MIN_SAFE_DISTANCE = 0.01f;  // nm (1% of typical sigma)
const float MAX_SAFE_ENERGY = 1e6f;     // kJ/mol

/**
 * @brief Calculate LJ and Coulomb energy with safety checks
 * 
 * !!! IMPORTANT: Zero distance handling strategy !!!
 * 1. For r < MIN_SAFE_DISTANCE:
 *    - Replace actual distance with MIN_SAFE_DISTANCE
 *    - Provides continuous potential without singularity
 *    - Allows MC moves through high-energy regions
 * 
 * 2. Energy capping:
 *    - Limits maximum repulsion to MAX_SAFE_ENERGY
 *    - Prevents exp(−βE) underflow in Metropolis
 *    - Maintains numerical stability of MC sampling
 * 
 * This approach:
 * - !!! Avoids infinite energies at r = 0
 * - !!! Keeps energy continuous and differentiable
 * - !!! Allows MC sampling of close contacts
 * - !!! Prevents numerical instabilities in simulation
 */
inline std::pair<float, float> calcPairEnergy(float r2, float sigma, float eps, float q1, float q2) {
    // !!! CRITICAL: Apply soft core potential for very small distances
    // This prevents infinities and numerical instabilities
    if (r2 < MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE) {
        r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;  // !!! Replace with safe minimum
    }
    
    float r = std::sqrt(r2);  // nm
    
    // Calculate LJ energy: V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    float sigma_r = sigma / r;
    float term6 = std::pow(sigma_r, 6);
    float term12 = term6 * term6;
    float vdw_energy = 4.0f * eps * (term12 - term6);  // kJ/mol
    
    // Calculate Coulomb energy: V_C = k_c * q1*q2/r
    float elec_energy = COULOMB * q1 * q2 / r;  // kJ/mol
    
    // !!! CRITICAL: Apply energy capping for numerical stability
    // First cap individual terms
    vdw_energy = std::min(vdw_energy, MAX_SAFE_ENERGY);
    vdw_energy = std::max(vdw_energy, -MAX_SAFE_ENERGY);
    elec_energy = std::min(elec_energy, MAX_SAFE_ENERGY);
    elec_energy = std::max(elec_energy, -MAX_SAFE_ENERGY);
    
    // !!! CRITICAL: Also cap total energy
    float total_energy = vdw_energy + elec_energy;
    if (total_energy > MAX_SAFE_ENERGY) {
        // Scale both components proportionally
        float scale = MAX_SAFE_ENERGY / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
    } else if (total_energy < -MAX_SAFE_ENERGY) {
        // Scale both components proportionally
        float scale = -MAX_SAFE_ENERGY / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
    }
    
    return {vdw_energy, elec_energy};
}

/**
 * @brief Compute non-bonded energies for movement residues without PBC
 * 
 * Calculates Lennard-Jones and Coulomb interactions between:
 * 1. Movement residues and all other active residues
 * 2. Each atom pair between residues
 * 
 * Energy components per residue:
 * 1. Lennard-Jones (kJ/mol):
 *    V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
 *    - ε: well depth (kJ/mol)
 *    - σ: distance at zero energy (nm)
 *    - r: inter-atomic distance (nm)
 * 
 * 2. Coulomb (kJ/mol):
 *    V_C = k_c * (q₁q₂/r)
 *    - k_c: Coulomb constant (138.935458 kJ·nm/mol/e²)
 *    - q₁,q₂: atomic charges (e)
 *    - r: inter-atomic distance (nm)
 * 
 * Uses soft core potential and energy capping for numerical stability.
 */
void computeNaiveNonbondedEnergy(model::MCState& state) {
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;

    // Validate force field parameter array sizes
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
        // Process only active movement residues
        for (int i = movementInfo.startIndex; 
             i < movementInfo.startIndex + movementInfo.activeCount; ++i) {
            if (!residues[i].active) continue;
            
            // For each atom in the movement residue
            for (int atom_i = residues[i].atomStart; 
                 atom_i < residues[i].atomStart + residues[i].atomCount; 
                 ++atom_i) {
                int moveType = atoms[atom_i].type;
                
                // Find movement type index in movementAtomTypes array
                int mi = -1;
                for (int k = 0; k < state.numMovementAtomTypes; ++k) {
                    if (state.movementAtomTypes[k] == moveType) {
                        mi = k;
                        break;
                    }
                }
                
                // Validate movement atom type
                if (mi < 0) {
                    std::stringstream ss;
                    ss << "Movement atom type " << moveType << " not found in movementAtomTypes "
                       << "for atom " << atom_i << " in residue " << i 
                       << " (" << movementInfo.resName << ")";
                    throw std::runtime_error(ss.str());
                }
                
                // Validate movement type index
                if (mi >= forcefield.numMovementTypes) {
                    std::stringstream ss;
                    ss << "Movement type index " << mi << " out of range. "
                       << "Maximum allowed index is " << (forcefield.numMovementTypes - 1)
                       << " for atom " << atom_i << " in residue " << i;
                    throw std::runtime_error(ss.str());
                }
                
                // Compute interaction with atoms in all other active residues
                for (int j = 0; j < state.activeResidueCount; ++j) {
                    // Skip if:
                    // 1. Same residue (avoid self-interaction)
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
                        
                        // Calculate distance between atoms (nm)
                        float dx = atoms[atom_j].x - atoms[atom_i].x;  // nm
                        float dy = atoms[atom_j].y - atoms[atom_i].y;  // nm
                        float dz = atoms[atom_j].z - atoms[atom_i].z;  // nm
                        float r2 = dx*dx + dy*dy + dz*dz;  // nm²
                        
                        // Get force field parameters
                        int param_index = mi * forcefield.numTotalTypes + resType;
                        float eps = forcefield.ljEps[param_index];     // kJ/mol
                        float sigma = forcefield.ljSigma[param_index]; // nm
                        float q1 = atoms[atom_i].charge;  // e
                        float q2 = atoms[atom_j].charge;  // e
                        
                        // Calculate energies with safety checks
                        auto [vdw_energy, elec_energy] = calcPairEnergy(r2, sigma, eps, q1, q2);
                        
                        // Add energies to the movement residue
                        residues[i].energy_vdw += vdw_energy;    // kJ/mol
                        residues[i].energy_elec += elec_energy;  // kJ/mol
                        
                        if (debug_output) {
                            float r = std::sqrt(r2);
                            platform::log(LogLevel::DEBUG,
                                "Interaction between residues ", i, "(", movementInfo.resName, ") and ", j, "\n",
                                "  Atoms: ", atom_i, "(type ", moveType, ") - ", 
                                atom_j, "(type ", resType, ")\n",
                                "  Distance: ", r, " nm\n",
                                "  Parameters: eps=", eps, " kJ/mol, sigma=", sigma, " nm\n",
                                "  Energies: vdw=", vdw_energy, " elec=", elec_energy, " kJ/mol\n",
                                "  Cumulative residue ", i, " energies: vdw=", residues[i].energy_vdw,
                                " elec=", residues[i].energy_elec, " kJ/mol");
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

