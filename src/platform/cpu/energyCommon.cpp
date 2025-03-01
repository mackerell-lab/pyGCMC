// src/platform/cpu/energyCommon.cpp

#include "energyCommon.hpp"
#include "energyDirect.hpp"
#include "energyEwald.hpp"
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
const float COULOMB = 138.935456f;
const float MIN_SAFE_DISTANCE = 0.01f;  // nm (1% of typical sigma)
const float MAX_SAFE_ENERGY = 1e6f;     // kJ/mol

// 调试标志
bool energy_debug_output = false;

/**
 * @brief 统一系统能量计算接口
 * 
 * @param state MC状态
 * @param method 能量计算方法（DIRECT或EWALD）
 * @param use_cutoff 是否使用截断
 * @param use_pbc 是否使用周期性边界条件
 */
void computeSystemEnergy(model::MCState& state, 
                         EnergyMethod method,
                         bool use_cutoff, 
                         bool use_pbc) {
    // 验证周期性边界条件的必要参数
    if (use_pbc) {
        validateBox(state.info.box, use_cutoff ? state.info.cutoff : 0.0f);
    }
    
    // 根据计算方法选择不同的实现
    switch (method) {
        case EnergyMethod::DIRECT:
            // 使用直接计算方法
            computeSystemEnergyDirect(state, use_cutoff, use_pbc);
            break;
            
        case EnergyMethod::EWALD:
            // 使用Ewald方法（需要周期性边界条件）
            if (!use_pbc) {
                throw std::runtime_error("Ewald method requires periodic boundary conditions");
            }
            computeSystemEnergyEwald(state);
            break;
    }
}

/**
 * @brief 统一运动残基能量计算接口
 * 
 * @param state MC状态
 * @param method 能量计算方法（DIRECT或EWALD）
 * @param use_cutoff 是否使用截断
 * @param use_pbc 是否使用周期性边界条件
 */
void computeMovementEnergy(model::MCState& state, 
                          EnergyMethod method,
                          bool use_cutoff, 
                          bool use_pbc) {
    // 验证周期性边界条件的必要参数
    if (use_pbc) {
        validateBox(state.info.box, use_cutoff ? state.info.cutoff : 0.0f);
    }
    
    // 根据计算方法选择不同的实现
    switch (method) {
        case EnergyMethod::DIRECT:
            // 使用直接计算方法
            computeMovementEnergyDirect(state, use_cutoff, use_pbc);
            break;
            
        case EnergyMethod::EWALD:
            // 使用Ewald方法（需要周期性边界条件）
            if (!use_pbc) {
                throw std::runtime_error("Ewald method requires periodic boundary conditions");
            }
            computeMovementEnergyEwald(state);
            break;
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 