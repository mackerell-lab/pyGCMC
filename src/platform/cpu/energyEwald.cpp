// src/platform/cpu/energyEwald.cpp

#include "energyEwald.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>  // For output formatting
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

// Define global variables
EwaldParams ewald_params;
const float TWO_PI = 2.0f * M_PI;

void EwaldParams::initializeTables(double cutoff) {
    this->cutoff = cutoff;
    ewaldDX = cutoff/NUM_TABLE_POINTS;
    ewaldDXInv = 1.0/ewaldDX;
    erfcDXInv = 1.0/(ewaldDX*alpha);
    
    erfcTable.resize(NUM_TABLE_POINTS + 4);
    ewaldScaleTable.resize(NUM_TABLE_POINTS + 4);
    
    for(int i = 0; i < NUM_TABLE_POINTS + 4; i++) {
        double r = i * ewaldDX;
        double alphaR = alpha * r;
        erfcTable[i] = std::erfc(alphaR);
        ewaldScaleTable[i] = erfcTable[i] + TWO_OVER_SQRT_PI * alphaR * std::exp(-alphaR * alphaR);
    }
}

void EwaldParams::initializeExpIkrTable(int numAtoms) {
    maxK = std::max(kmax[0], std::max(kmax[1], kmax[2]));
    expIkrTable.resize(maxK * numAtoms * 3);
    expIkrXY.resize(numAtoms);
}

double EwaldParams::erfcApprox(double r) const {
    double x = r * erfcDXInv;
    int index = std::min(static_cast<int>(x), NUM_TABLE_POINTS);
    double coeff2 = x - index;
    double coeff1 = 1.0 - coeff2;
    return coeff1 * erfcTable[index] + coeff2 * erfcTable[index + 1];
}

double EwaldParams::ewaldScaleApprox(double r) const {
    double x = r * ewaldDXInv;
    int index = std::min(static_cast<int>(x), NUM_TABLE_POINTS);
    double coeff2 = x - index;
    double coeff1 = 1.0 - coeff2;
    return coeff1 * ewaldScaleTable[index] + coeff2 * ewaldScaleTable[index + 1];
}

void autoAdjustParameters(double error_tolerance, double cutoff_distance, const double box[3]) {
    // 检查cutoff是否小于盒子长度的一半
    double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
    if (cutoff_distance >= 0.5 * minBoxSize) {
        throw std::runtime_error("Cutoff distance must be less than half the smallest box dimension");
    }
    
    // Calculate optimal alpha based on error tolerance and cutoff
    ewald_params.alpha = std::sqrt(-std::log(2.0 * error_tolerance)) / cutoff_distance;
    
    // Calculate optimal kmax for each dimension
    double kmax_float = 2.0 * ewald_params.alpha * minBoxSize * 
                      std::sqrt(-std::log(2.0 * error_tolerance));
    
    for(int i = 0; i < 3; i++) {
        ewald_params.kmax[i] = static_cast<int>(std::ceil(kmax_float * minBoxSize/box[i]));
    }
    
    // Initialize lookup tables
    ewald_params.initializeTables(cutoff_distance);
    ewald_params.initialized = true;
}

/**
 * @brief 设置Ewald计算参数
 * 
 * @param alpha Ewald分离参数 (nm^-1)
 * @param kmax Maximum reciprocal space wave vectors
 * @param tolerance Precision control
 */
void setEwaldParameters(double alpha, const int kmax[3], double tolerance) {
    ewald_params.alpha = alpha;
    for(int i = 0; i < 3; i++) {
        ewald_params.kmax[i] = kmax[i];
    }
    ewald_params.tolerance = tolerance;
    ewald_params.initialized = true;
}

/**
 * @brief 计算Ewald实空间部分的能量
 * 
 * The real space part includes:
 * 1. van der Waals interactions (same as direct calculation)
 * 2. short-range Coulomb interactions (erfc(αr)/r)
 */
inline std::pair<double, double> calcPairEnergyEwald(
    double r2, double sigma, double eps, double q1, double q2) {
    
    // 应用最小安全距离
    if (r2 < MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE) {
        r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;
    }
    
    double r = std::sqrt(r2);
    
    // VDW能量计算保持不变
    double sigma_r2 = (sigma * sigma) / r2;
    double sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;
    double sigma_r12 = sigma_r6 * sigma_r6;
    double vdw_energy = 4.0 * eps * (sigma_r12 - sigma_r6);
    
    // Remove 0.5 factor since outer loop already handles i<j pairs
    double erfc_term = ewald_params.erfcApprox(r);
    double elec_energy = COULOMB * q1 * q2 * erfc_term / r;
    
    // 应用能量上限
    const double max_safe_energy = static_cast<double>(MAX_SAFE_ENERGY);
    vdw_energy = std::min(std::max(vdw_energy, -max_safe_energy), max_safe_energy);
    elec_energy = std::min(std::max(elec_energy, -max_safe_energy), max_safe_energy);
    
    return {vdw_energy, elec_energy};
}

/**
 * @brief 计算Ewald倒空间部分的能量
 * 
 * @param state 系统状态
 * @param movement_only 是否只计算movement residues
 * @return 倒空间总能量
 */
double computeReciprocalEnergy(model::MCState& state, bool movement_only) {
    const auto& box = state.info.box;
    const auto& atoms = state.atoms;
    double volume = box[0] * box[1] * box[2];
    
    // Calculate total charge
    double totalCharge = 0.0;
    for(const auto& atom : atoms) {
        totalCharge += static_cast<double>(atom.charge);
    }
    
    typedef std::complex<double> Complex;
    const double TWO_PI = 2.0 * M_PI;
    // Standard reciprocal coefficient (2π/V)
    const double recipCoeff = COULOMB * TWO_PI / volume;
    const double factorEwald = -1.0 / (4.0 * ewald_params.alpha * ewald_params.alpha);
    
    double total_energy = 0.0;
    
    // Standard k-space summation over all k-vectors
    for (int rx = -ewald_params.kmax[0]; rx <= ewald_params.kmax[0]; rx++) {
        double kx = rx * TWO_PI / box[0];
        
        for (int ry = -ewald_params.kmax[1]; ry <= ewald_params.kmax[1]; ry++) {
            double ky = ry * TWO_PI / box[1];
            
            for (int rz = -ewald_params.kmax[2]; rz <= ewald_params.kmax[2]; rz++) {
                // Skip k = 0
                if (rx == 0 && ry == 0 && rz == 0) continue;
                
                double kz = rz * TWO_PI / box[2];
                double k2 = kx*kx + ky*ky + kz*kz;
                
                Complex structureFactor(0.0, 0.0);
                for (int n = 0; n < static_cast<int>(atoms.size()); n++) {
                    if (movement_only) {
                        bool in_movement = false;
                        for (const auto& movementInfo : state.movementResidues) {
                            if (n >= movementInfo.startIndex && 
                                n < movementInfo.startIndex + movementInfo.activeCount) {
                                in_movement = true;
                                break;
                            }
                        }
                        if (!in_movement) continue;
                    }
                    
                    double kdotr = kx*static_cast<double>(atoms[n].x) + 
                                 ky*static_cast<double>(atoms[n].y) + 
                                 kz*static_cast<double>(atoms[n].z);
                    Complex phase(std::cos(kdotr), std::sin(kdotr));
                    structureFactor += static_cast<double>(atoms[n].charge) * phase;
                }
                
                double ak = std::exp(k2 * factorEwald) / k2;
                double structureFactorNorm = std::norm(structureFactor);
                
                // No symmetry factor needed - we sum over all k-vectors
                total_energy += recipCoeff * ak * structureFactorNorm;
            }
        }
    }
    
    return total_energy;
}

/**
 * @brief 计算自能校正项
 */
double computeSelfEnergy(model::MCState& state, bool movement_only) {
    double self_energy = 0.0;
    double totalCharge = 0.0;
    
    if(movement_only) {
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex; 
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(!state.residues[i].active) continue;
                
                for(int j = state.residues[i].atomStart;
                    j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                    double charge = state.atoms[j].charge;
                    // Standard self-energy term
                    self_energy -= charge * charge;
                    totalCharge += charge;
                }
            }
        }
    } else {
        for(int r = 0; r < state.activeResidueCount; r++) {
            if(!state.residues[r].active) continue;
            
            for(int i = state.residues[r].atomStart;
                i < state.residues[r].atomStart + state.residues[r].atomCount; i++) {
                double charge = state.atoms[i].charge;
                // Standard self-energy term
                self_energy -= charge * charge;
                totalCharge += charge;
            }
        }
    }
    
    // Standard self-energy prefactor
    double baseEnergy = self_energy * COULOMB * ewald_params.alpha / std::sqrt(M_PI);
    
    // Background correction for non-zero net charge
    if (!movement_only && std::abs(totalCharge) > 1e-10) {
        double volume = state.info.box[0] * state.info.box[1] * state.info.box[2];
        double backgroundCorrection = -COULOMB * M_PI * totalCharge * totalCharge / 
            (2.0 * volume * ewald_params.alpha * ewald_params.alpha);
        baseEnergy += backgroundCorrection;
    }
    
    return baseEnergy;
}

void checkSystemNeutrality(const model::MCState& state) {
    float totalCharge = 0.0f;
    for(const auto& atom : state.atoms) {
        totalCharge += atom.charge;
    }
    if(std::abs(totalCharge) > 1e-6f) {
        throw std::runtime_error("Ewald summation requires neutral system");
    }
}

/**
 * @brief 使用Ewald方法计算系统能量
 */
void computeSystemEnergyEwald(model::MCState& state) {
    if (!ewald_params.initialized) {
        throw std::runtime_error("Ewald parameters not initialized");
    }
    
    // 检查PBC和cutoff条件
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Ewald method requires periodic boundary conditions");
    }
    
    float minBoxSize = std::min(state.info.box[0], std::min(state.info.box[1], state.info.box[2]));
    if (ewald_params.cutoff >= 0.5f * minBoxSize) {
        throw std::runtime_error("Cutoff distance must be less than half the smallest box dimension");
    }
    
    // 实空间部分
    computeSystemEnergyCutoff(state);
    
    // 倒空间部分
    double recip_energy = computeReciprocalEnergy(state, false);
    
    // 自能校正
    double self_energy = computeSelfEnergy(state, false);
    
    // 分配长程能量到residues
    int active_count = 0;
    for(const auto& residue : state.residues) {
        if(residue.active) active_count++;
    }
    
    if(active_count > 0) {
        double energy_per_residue = (recip_energy + self_energy) / active_count;
        for(auto& residue : state.residues) {
            if(residue.active) {
                residue.energy_elec += energy_per_residue;
            }
        }
    }
}

/**
 * @brief 使用Ewald方法计算movement residues的能量
 */
void computeMovementEnergyEwald(model::MCState& state) {
    if (!ewald_params.initialized) {
        throw std::runtime_error("Ewald parameters not initialized");
    }
    
    // 检查PBC和cutoff条件
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Ewald method requires periodic boundary conditions");
    }
    
    float minBoxSize = std::min({state.info.box[0], state.info.box[1], state.info.box[2]});
    if (ewald_params.cutoff >= 0.5f * minBoxSize) {
        throw std::runtime_error("Cutoff distance must be less than half the smallest box dimension");
    }
    
    // 实空间部分
    computeMovementEnergyCutoff(state);
    
    // 倒空间部分
    double recip_energy = computeReciprocalEnergy(state, true);
    
    // 自能校正
    double self_energy = computeSelfEnergy(state, true);
    
    // 分配长程能量到movement residues
    int movement_count = 0;
    for(const auto& movementInfo : state.movementResidues) {
        movement_count += movementInfo.activeCount;
    }
    
    if(movement_count > 0) {
        double energy_per_residue = (recip_energy + self_energy) / movement_count;
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex;
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(state.residues[i].active) {
                    state.residues[i].energy_elec += energy_per_residue;
                }
            }
        }
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
