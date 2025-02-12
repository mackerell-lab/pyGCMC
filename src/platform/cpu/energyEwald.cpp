// src/platform/cpu/energyEwald.cpp
#include "energy.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace pygcmc {
namespace platform {
namespace cpu {

// 定义全局变量
EwaldParams ewald_params;

/**
 * @brief 设置Ewald计算参数
 * 
 * @param alpha Ewald分离参数 (nm^-1)
 * @param kmax 倒空间最大波矢
 * @param tolerance 精度控制
 */
void setEwaldParameters(float alpha, const int kmax[3], float tolerance) {
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
 * 实空间部分包含:
 * 1. 范德华相互作用 (与直接计算相同)
 * 2. 短程库仑相互作用 (erfc(αr)/r)
 */
inline std::pair<float, float> calcPairEnergyEwald(
    float r2, float sigma, float eps, float q1, float q2) {
    
    // 应用最小安全距离
    if (r2 < MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE) {
        r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;
    }
    
    float r = std::sqrt(r2);
    
    // VDW能量计算保持不变
    float sigma_r2 = (sigma * sigma) / r2;
    float sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;
    float sigma_r12 = sigma_r6 * sigma_r6;
    float vdw_energy = 4.0f * eps * (sigma_r12 - sigma_r6);
    
    // Ewald实空间静电能
    float alpha_r = ewald_params.alpha * r;
    float erfc_term = std::erfc(alpha_r);
    float elec_energy = COULOMB * q1 * q2 * erfc_term / r;
    
    // 应用能量上限
    vdw_energy = std::min(std::max(vdw_energy, -MAX_SAFE_ENERGY), MAX_SAFE_ENERGY);
    elec_energy = std::min(std::max(elec_energy, -MAX_SAFE_ENERGY), MAX_SAFE_ENERGY);
    
    return {vdw_energy, elec_energy};
}

/**
 * @brief 计算Ewald倒空间部分的能量
 * 
 * @param state 系统状态
 * @param movement_only 是否只计算movement residues
 * @return 倒空间总能量
 */
float computeReciprocalEnergy(model::MCState& state, bool movement_only) {
    const auto& box = state.info.box;
    const auto& atoms = state.atoms;
    float volume = box[0] * box[1] * box[2];
    
    // 计算倒空间基矢
    float recip_vec[3][3];
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            recip_vec[i][j] = (i == j) ? 2.0f * M_PI / box[i] : 0.0f;
        }
    }
    
    float total_energy = 0.0f;
    
    // 倒空间求和
    for(int kx = -ewald_params.kmax[0]; kx <= ewald_params.kmax[0]; kx++) {
        for(int ky = -ewald_params.kmax[1]; ky <= ewald_params.kmax[1]; ky++) {
            for(int kz = -ewald_params.kmax[2]; kz <= ewald_params.kmax[2]; kz++) {
                if(kx == 0 && ky == 0 && kz == 0) continue;
                
                float k[3] = {
                    kx * recip_vec[0][0],
                    ky * recip_vec[1][1],
                    kz * recip_vec[2][2]
                };
                
                float k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
                float k_factor = 2.0f * M_PI / volume * 
                               std::exp(-k2/(4.0f*ewald_params.alpha*ewald_params.alpha)) / k2;
                
                // 计算结构因子
                float struct_real = 0.0f, struct_imag = 0.0f;
                
                if(movement_only) {
                    // 只计算movement residues的贡献
                    for(const auto& movementInfo : state.movementResidues) {
                        for(int i = movementInfo.startIndex; 
                            i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                            if(!state.residues[i].active) continue;
                            
                            for(int j = state.residues[i].atomStart;
                                j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                                float kdotr = k[0]*atoms[j].x + k[1]*atoms[j].y + k[2]*atoms[j].z;
                                float charge = atoms[j].charge;
                                struct_real += charge * std::cos(kdotr);
                                struct_imag += charge * std::sin(kdotr);
                            }
                        }
                    }
                } else {
                    // 计算所有active residues的贡献
                    for(int r = 0; r < state.activeResidueCount; r++) {
                        if(!state.residues[r].active) continue;
                        
                        for(int i = state.residues[r].atomStart;
                            i < state.residues[r].atomStart + state.residues[r].atomCount; i++) {
                            float kdotr = k[0]*atoms[i].x + k[1]*atoms[i].y + k[2]*atoms[i].z;
                            float charge = atoms[i].charge;
                            struct_real += charge * std::cos(kdotr);
                            struct_imag += charge * std::sin(kdotr);
                        }
                    }
                }
                
                total_energy += k_factor * (struct_real*struct_real + struct_imag*struct_imag);
            }
        }
    }
    
    return total_energy * COULOMB;
}

/**
 * @brief 计算自能校正项
 */
float computeSelfEnergy(model::MCState& state, bool movement_only) {
    float self_energy = 0.0f;
    
    if(movement_only) {
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex; 
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(!state.residues[i].active) continue;
                
                for(int j = state.residues[i].atomStart;
                    j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                    float charge = state.atoms[j].charge;
                    self_energy -= charge * charge;
                }
            }
        }
    } else {
        for(int r = 0; r < state.activeResidueCount; r++) {
            if(!state.residues[r].active) continue;
            
            for(int i = state.residues[r].atomStart;
                i < state.residues[r].atomStart + state.residues[r].atomCount; i++) {
                float charge = state.atoms[i].charge;
                self_energy -= charge * charge;
            }
        }
    }
    
    return self_energy * COULOMB * ewald_params.alpha / std::sqrt(M_PI);
}

/**
 * @brief 使用Ewald方法计算系统能量
 */
void computeSystemEnergyEwald(model::MCState& state) {
    if (!ewald_params.initialized) {
        throw std::runtime_error("Ewald parameters not initialized");
    }
    
    // 检查PBC
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Ewald method requires periodic boundary conditions");
    }
    
    // 实空间部分 (使用现有的cutoff PBC计算)
    computeSystemEnergyCutoff(state);
    
    // 倒空间部分
    float recip_energy = computeReciprocalEnergy(state, false);
    
    // 自能校正
    float self_energy = computeSelfEnergy(state, false);
    
    // 分配长程能量到residues
    int active_count = 0;
    for(const auto& residue : state.residues) {
        if(residue.active) active_count++;
    }
    
    if(active_count > 0) {
        float energy_per_residue = (recip_energy + self_energy) / active_count;
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
    
    // 检查PBC
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Ewald method requires periodic boundary conditions");
    }
    
    // 实空间部分
    computeMovementEnergyCutoff(state);
    
    // 倒空间部分
    float recip_energy = computeReciprocalEnergy(state, true);
    
    // 自能校正
    float self_energy = computeSelfEnergy(state, true);
    
    // 分配长程能量到movement residues
    int movement_count = 0;
    for(const auto& movementInfo : state.movementResidues) {
        movement_count += movementInfo.activeCount;
    }
    
    if(movement_count > 0) {
        float energy_per_residue = (recip_energy + self_energy) / movement_count;
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
