// src/platform/cpu/energy.cpp
#include "energy.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>  // 添加用于格式化输出

namespace pygcmc {
namespace platform {
namespace cpu {

// Debug flag to control output
static bool debug_output = true;  // 默认开启调试输出

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
    if (debug_output) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\n=== calcPairEnergy called ===";
        ss << "\nInput parameters:"
           << "\n  Distance² = " << r2 << " nm²"
           << "\n  Sigma = " << sigma << " nm"
           << "\n  Epsilon = " << eps << " kJ/mol"
           << "\n  q1 = " << q1 << " e"
           << "\n  q2 = " << q2 << " e";
        platform::log(LogLevel::DEBUG, ss.str());
    }

    if (r2 < MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE) {
        if (debug_output) {
            platform::log(LogLevel::DEBUG, "Distance below MIN_SAFE_DISTANCE, using r2 = ", 
                         MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE);
        }
        r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;
    }
    
    float r = std::sqrt(r2);
    
    if (debug_output) {
        platform::log(LogLevel::DEBUG, "Distance r = ", r, " nm");
    }
    
    // Calculate LJ energy: V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    float sigma_r = sigma / r;
    float term6 = std::pow(sigma_r, 6);
    float term12 = term6 * term6;
    float vdw_energy = 4.0f * eps * (term12 - term6);  // kJ/mol
    
    // Calculate Coulomb energy: V_C = k_c * q1*q2/r
    float elec_energy = COULOMB * q1 * q2 / r;  // kJ/mol

    if (debug_output) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\nEnergy calculation details:";
        ss << "\n  sigma/r = " << sigma_r;
        ss << "\n  (sigma/r)^6 = " << term6;
        ss << "\n  (sigma/r)^12 = " << term12;
        ss << "\n  4*epsilon = " << (4.0f * eps);
        ss << "\n  VDW term = " << (term12 - term6);
        ss << "\n  COULOMB constant = " << COULOMB;
        ss << "\n  q1*q2 = " << (q1 * q2);
        ss << "\nInitial energies:";
        ss << "\n  VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Electrostatic energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // !!! CRITICAL: Apply energy capping for numerical stability
    // First cap individual terms
    vdw_energy = std::min(vdw_energy, MAX_SAFE_ENERGY);
    vdw_energy = std::max(vdw_energy, -MAX_SAFE_ENERGY);
    elec_energy = std::min(elec_energy, MAX_SAFE_ENERGY);
    elec_energy = std::max(elec_energy, -MAX_SAFE_ENERGY);
    
    // !!! CRITICAL: Also cap total energy
    float total_energy = vdw_energy + elec_energy;
    if (total_energy > MAX_SAFE_ENERGY) {
        float scale = MAX_SAFE_ENERGY / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (debug_output) {
            platform::log(LogLevel::DEBUG, "Total energy exceeded MAX_SAFE_ENERGY, scaled by ", scale);
        }
    } else if (total_energy < -MAX_SAFE_ENERGY) {
        float scale = -MAX_SAFE_ENERGY / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (debug_output) {
            platform::log(LogLevel::DEBUG, "Total energy below -MAX_SAFE_ENERGY, scaled by ", scale);
        }
    }

    if (debug_output) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\n=== calcPairEnergy returning ===";
        ss << "\n  Final VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Final Electrostatic energy = " << elec_energy << " kJ/mol";
        ss << "\n  Total energy = " << (vdw_energy + elec_energy) << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
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
    if (debug_output) {
        platform::log(LogLevel::DEBUG, "\n=== Starting nonbonded energy calculation ===");
        platform::log(LogLevel::DEBUG, "System state info:");
        platform::log(LogLevel::DEBUG, "  Active residue count: ", state.activeResidueCount);
        platform::log(LogLevel::DEBUG, "  Number of movement residues: ", state.movementResidues.size());
        platform::log(LogLevel::DEBUG, "  Number of movement atom types: ", state.numMovementAtomTypes);
        platform::log(LogLevel::DEBUG, "Force field info:");
        platform::log(LogLevel::DEBUG, "  Number of movement types: ", state.forcefield.numMovementTypes);
        platform::log(LogLevel::DEBUG, "  Number of total types: ", state.forcefield.numTotalTypes);
        platform::log(LogLevel::DEBUG, "  Size of LJ parameters: ", state.forcefield.ljEps.size());
    }

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
        if (debug_output) {
            platform::log(LogLevel::DEBUG, "Reset energies for residue at atom start ", residue.atomStart);
        }
    }

    // Iterate through all movement molecule groups
    for (const auto& movementInfo : state.movementResidues) {
        if (debug_output) {
            platform::log(LogLevel::DEBUG, "\nProcessing movement residue group: ", movementInfo.resName);
            platform::log(LogLevel::DEBUG, "  Start index: ", movementInfo.startIndex);
            platform::log(LogLevel::DEBUG, "  Active count: ", movementInfo.activeCount);
            platform::log(LogLevel::DEBUG, "  Total count: ", movementInfo.totalCount);
        }

        // Process only active movement residues
        for (int i = movementInfo.startIndex; 
             i < movementInfo.startIndex + movementInfo.activeCount; ++i) {
            if (!residues[i].active) continue;
            
            if (debug_output) {
                platform::log(LogLevel::DEBUG, "\nProcessing movement residue ", i);
                platform::log(LogLevel::DEBUG, "  Atom start: ", residues[i].atomStart);
                platform::log(LogLevel::DEBUG, "  Atom count: ", residues[i].atomCount);
            }

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
                
                if (debug_output) {
                    platform::log(LogLevel::DEBUG, "\nProcessing movement atom ", atom_i);
                    platform::log(LogLevel::DEBUG, "  Type: ", moveType);
                    platform::log(LogLevel::DEBUG, "  Movement type index: ", mi);
                    platform::log(LogLevel::DEBUG, "  Position: (", 
                                atoms[atom_i].x, ", ", 
                                atoms[atom_i].y, ", ", 
                                atoms[atom_i].z, ") nm");
                }
                
                // Compute interaction with atoms in all other active residues
                for (int j = 0; j < state.activeResidueCount; ++j) {
                    // Skip if:
                    // 1. Same residue (avoid self-interaction)
                    // 2. Residue is inactive
                    if (j == i || !residues[j].active) continue;
                    
                    if (debug_output) {
                        platform::log(LogLevel::DEBUG, "\n  Interacting with residue ", j);
                        platform::log(LogLevel::DEBUG, "    Atom start: ", residues[j].atomStart);
                        platform::log(LogLevel::DEBUG, "    Atom count: ", residues[j].atomCount);
                    }

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
                        
                        if (debug_output) {
                            platform::log(LogLevel::DEBUG, "\n    Interacting with atom ", atom_j);
                            platform::log(LogLevel::DEBUG, "      Type: ", resType);
                            platform::log(LogLevel::DEBUG, "      Position: (", 
                                        atoms[atom_j].x, ", ", 
                                        atoms[atom_j].y, ", ", 
                                        atoms[atom_j].z, ") nm");
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
                        
                        if (debug_output) {
                            std::stringstream ss;
                            ss << std::fixed << std::setprecision(6);
                            ss << "\n      Interaction parameters:";
                            ss << "\n        Parameter index: " << param_index;
                            ss << "\n        Epsilon: " << eps << " kJ/mol";
                            ss << "\n        Sigma: " << sigma << " nm";
                            ss << "\n        Charges: q1=" << q1 << "e, q2=" << q2 << "e";
                            ss << "\n        Distance vector: (" << dx << ", " << dy << ", " << dz << ") nm";
                            ss << "\n        Distance²: " << r2 << " nm²";
                            platform::log(LogLevel::DEBUG, ss.str());
                        }
                        
                        // Calculate energies with safety checks
                        auto [vdw_energy, elec_energy] = calcPairEnergy(r2, sigma, eps, q1, q2);
                        
                        // Add energies to the movement residue
                        residues[i].energy_vdw += vdw_energy;
                        residues[i].energy_elec += elec_energy;
                        residues[j].energy_vdw += vdw_energy;
                        residues[j].energy_elec += elec_energy;

                        if (debug_output) {
                            std::stringstream ss;
                            ss << std::fixed << std::setprecision(6);
                            ss << "\n      Cumulative residue " << i << " energies:";
                            ss << "\n        VDW: " << residues[i].energy_vdw << " kJ/mol";
                            ss << "\n        Electrostatic: " << residues[i].energy_elec << " kJ/mol";
                            ss << "\n        Total: " << (residues[i].energy_vdw + residues[i].energy_elec) << " kJ/mol";
                            platform::log(LogLevel::DEBUG, ss.str());
                        }
                    }
                }
            }
        }
    }

    if (debug_output) {
        platform::log(LogLevel::DEBUG, "\n=== Final energies for all residues ===");
        for (size_t i = 0; i < residues.size(); ++i) {
            if (residues[i].active) {
                platform::log(LogLevel::DEBUG, "Residue ", i, ":");
                platform::log(LogLevel::DEBUG, "  VDW energy: ", residues[i].energy_vdw, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Electrostatic energy: ", residues[i].energy_elec, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Total energy: ", 
                            (residues[i].energy_vdw + residues[i].energy_elec), " kJ/mol");
            }
        }
        platform::log(LogLevel::DEBUG, "\n=== Completed nonbonded energy calculation ===");
    }
}

void computeAllNonbondedEnergy(model::MCState& state) {
    if (debug_output) {
        platform::log(LogLevel::DEBUG, "\n=== Starting all residues nonbonded energy calculation ===");
        platform::log(LogLevel::DEBUG, "System state info:");
        platform::log(LogLevel::DEBUG, "  Active residue count: ", state.activeResidueCount);
    }

    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;

    // 验证力场参数数组大小
    size_t expected_size = static_cast<size_t>(forcefield.numTotalTypes) * 
                          static_cast<size_t>(forcefield.numTotalTypes);
    if (forcefield.ljEps.size() != expected_size || forcefield.ljSigma.size() != expected_size) {
        std::stringstream ss;
        ss << "Force field parameters array size mismatch. Expected size "
           << expected_size
           << " (numTotalTypes * numTotalTypes), but got eps=" << forcefield.ljEps.size()
           << " sigma=" << forcefield.ljSigma.size();
        throw std::runtime_error(ss.str());
    }

    // 重置所有residue的能量
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        residue.energy_elec = 0.0f;
    }

    // 遍历所有active residues
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (!residues[i].active) continue;
        
        // 遍历当前residue的所有原子
        for (int atom_i = residues[i].atomStart; 
             atom_i < residues[i].atomStart + residues[i].atomCount; 
             ++atom_i) {
            int type_i = atoms[atom_i].type;
            
            // 验证原子类型
            if (type_i >= forcefield.numTotalTypes) {
                std::stringstream ss;
                ss << "Atom type " << type_i << " out of range. "
                   << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
                   << " for atom " << atom_i << " in residue " << i;
                throw std::runtime_error(ss.str());
            }

            // 与其他active residues的原子计算相互作用
            for (int j = i + 1; j < state.activeResidueCount; ++j) {
                if (!residues[j].active) continue;

                // 遍历另一个residue的所有原子
                for (int atom_j = residues[j].atomStart;
                     atom_j < residues[j].atomStart + residues[j].atomCount;
                     ++atom_j) {
                    int type_j = atoms[atom_j].type;
                    
                    // 验证原子类型
                    if (type_j >= forcefield.numTotalTypes) {
                        std::stringstream ss;
                        ss << "Atom type " << type_j << " out of range. "
                           << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
                           << " for atom " << atom_j << " in residue " << j;
                        throw std::runtime_error(ss.str());
                    }

                    // 计算原子间距离
                    float dx = atoms[atom_j].x - atoms[atom_i].x;
                    float dy = atoms[atom_j].y - atoms[atom_i].y;
                    float dz = atoms[atom_j].z - atoms[atom_i].z;
                    float r2 = dx*dx + dy*dy + dz*dz;

                    // 获取力场参数
                    int param_index = type_i * forcefield.numTotalTypes + type_j;
                    float eps = forcefield.ljEps[param_index];
                    float sigma = forcefield.ljSigma[param_index];
                    float q1 = atoms[atom_i].charge;
                    float q2 = atoms[atom_j].charge;

                    // 计算能量
                    auto [vdw_energy, elec_energy] = calcPairEnergy(r2, sigma, eps, q1, q2);

                    // 将能量完整地加到每个residue上
                    residues[i].energy_vdw += vdw_energy;
                    residues[i].energy_elec += elec_energy;
                    residues[j].energy_vdw += vdw_energy;
                    residues[j].energy_elec += elec_energy;
                }
            }
        }
    }

    if (debug_output) {
        platform::log(LogLevel::DEBUG, "\n=== Final energies for all residues ===");
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (residues[i].active) {
                platform::log(LogLevel::DEBUG, "Residue ", i, ":");
                platform::log(LogLevel::DEBUG, "  VDW energy: ", residues[i].energy_vdw, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Electrostatic energy: ", residues[i].energy_elec, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Total energy: ", 
                            (residues[i].energy_vdw + residues[i].energy_elec), " kJ/mol");
            }
        }
        platform::log(LogLevel::DEBUG, "\n=== Completed all residues nonbonded energy calculation ===");
    }
}

void computeCutoffNonPeriodicEnergy(model::MCState& state) {
    if (debug_output) {
        platform::log(LogLevel::DEBUG, "\n=== Starting cutoff nonbonded energy calculation ===");
        platform::log(LogLevel::DEBUG, "System state info:");
        platform::log(LogLevel::DEBUG, "  Active residue count: ", state.activeResidueCount);
        platform::log(LogLevel::DEBUG, "  Cutoff distance: ", state.info.cutoff, " nm");
    }

    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;
    const float cutoff2 = state.info.cutoff * state.info.cutoff;  // 截断距离的平方

    // 验证力场参数数组大小
    size_t expected_size = static_cast<size_t>(forcefield.numTotalTypes) * 
                          static_cast<size_t>(forcefield.numTotalTypes);
    if (forcefield.ljEps.size() != expected_size || forcefield.ljSigma.size() != expected_size) {
        std::stringstream ss;
        ss << "Force field parameters array size mismatch. Expected size "
           << expected_size
           << " (numTotalTypes * numTotalTypes), but got eps=" << forcefield.ljEps.size()
           << " sigma=" << forcefield.ljSigma.size();
        throw std::runtime_error(ss.str());
    }

    // 重置所有residue的能量
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        residue.energy_elec = 0.0f;
    }

    // 遍历所有active residues
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (!residues[i].active) continue;
        
        // 遍历当前residue的所有原子
        for (int atom_i = residues[i].atomStart; 
             atom_i < residues[i].atomStart + residues[i].atomCount; 
             ++atom_i) {
            int type_i = atoms[atom_i].type;
            
            // 验证原子类型
            if (type_i >= forcefield.numTotalTypes) {
                std::stringstream ss;
                ss << "Atom type " << type_i << " out of range. "
                   << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
                   << " for atom " << atom_i << " in residue " << i;
                throw std::runtime_error(ss.str());
            }

            // 与其他active residues的原子计算相互作用
            for (int j = i + 1; j < state.activeResidueCount; ++j) {
                if (!residues[j].active) continue;

                // 遍历另一个residue的所有原子
                for (int atom_j = residues[j].atomStart;
                     atom_j < residues[j].atomStart + residues[j].atomCount;
                     ++atom_j) {
                    int type_j = atoms[atom_j].type;
                    
                    // 验证原子类型
                    if (type_j >= forcefield.numTotalTypes) {
                        std::stringstream ss;
                        ss << "Atom type " << type_j << " out of range. "
                           << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
                           << " for atom " << atom_j << " in residue " << j;
                        throw std::runtime_error(ss.str());
                    }

                    // 计算原子间距离
                    float dx = atoms[atom_j].x - atoms[atom_i].x;
                    float dy = atoms[atom_j].y - atoms[atom_i].y;
                    float dz = atoms[atom_j].z - atoms[atom_i].z;
                    float r2 = dx*dx + dy*dy + dz*dz;

                    // 检查是否在截断距离内
                    if (r2 > cutoff2) continue;

                    // 获取力场参数
                    int param_index = type_i * forcefield.numTotalTypes + type_j;
                    float eps = forcefield.ljEps[param_index];
                    float sigma = forcefield.ljSigma[param_index];
                    float q1 = atoms[atom_i].charge;
                    float q2 = atoms[atom_j].charge;

                    // 计算能量
                    auto [vdw_energy, elec_energy] = calcPairEnergy(r2, sigma, eps, q1, q2);

                    // 将能量完整地加到每个residue上
                    residues[i].energy_vdw += vdw_energy;
                    residues[i].energy_elec += elec_energy;
                    residues[j].energy_vdw += vdw_energy;
                    residues[j].energy_elec += elec_energy;
                }
            }
        }
    }

    if (debug_output) {
        platform::log(LogLevel::DEBUG, "\n=== Final energies for all residues ===");
        float total_vdw = 0.0f;
        float total_elec = 0.0f;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (residues[i].active) {
                total_vdw += residues[i].energy_vdw;
                total_elec += residues[i].energy_elec;
                platform::log(LogLevel::DEBUG, "Residue ", i, ":");
                platform::log(LogLevel::DEBUG, "  VDW energy: ", residues[i].energy_vdw, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Electrostatic energy: ", residues[i].energy_elec, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Total energy: ", 
                            (residues[i].energy_vdw + residues[i].energy_elec), " kJ/mol");
            }
        }
        platform::log(LogLevel::DEBUG, "\nTotal system energy (before /2):");
        platform::log(LogLevel::DEBUG, "  VDW: ", total_vdw, " kJ/mol");
        platform::log(LogLevel::DEBUG, "  Electrostatic: ", total_elec, " kJ/mol");
        platform::log(LogLevel::DEBUG, "  Total: ", (total_vdw + total_elec), " kJ/mol");
        platform::log(LogLevel::DEBUG, "\n=== Completed cutoff nonbonded energy calculation ===");
    }
}

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable) {
    debug_output = enable;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc

