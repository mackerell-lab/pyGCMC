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
 * @brief 配置 CHARMM 风格的平滑函数参数
 * 
 * @param state MC状态，将在其中更新switching参数
 * @param use_switching 是否启用平滑函数
 * @param r_on 内截断半径 (ctonnb)
 * @param r_off 外截断半径 (ctofnb)
 */
void setSwitchingFunction(model::MCState& state, bool use_switching, float r_on, float r_off) {
    // Input validation
    if (r_on >= r_off) {
        throw std::runtime_error("Invalid switching function parameters: r_on must be less than r_off");
    }
    if (r_on <= 0.0f || r_off <= 0.0f) {
        throw std::runtime_error("Invalid switching function parameters: radii must be positive");
    }

    // Set parameters in MCState
    state.info.use_switching = use_switching;
    state.info.r_on = r_on;
    state.info.r_off = r_off;
    
    if (energy_debug_output && use_switching) {
        platform::log(LogLevel::INFO, 
            "CHARMM switching function enabled: r_on=", r_on, " nm, r_off=", r_off, " nm");
    } else if (energy_debug_output) {
        platform::log(LogLevel::INFO, "CHARMM switching function disabled");
    }
}

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