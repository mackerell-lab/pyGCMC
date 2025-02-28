// src/platform/cpu/energy.hpp

#pragma once

// 首先包含通用接口
#include "energyCommon.hpp"

// 然后包含实现
#include "energyDirect.hpp"
#include "energyEwald.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief CPU平台能量计算模块
 * 
 * 本模块提供了CPU平台上的能量计算功能，包括：
 * 1. 直接计算方法（Direct）：使用显式求和计算范德华和库仑相互作用
 * 2. Ewald求和方法：针对周期性系统的长程静电相互作用进行优化计算
 * 
 * 基本用法示例：
 * 
 * // 直接计算（无截断、无PBC）：
 * computeSystemEnergy(state, EnergyMethod::DIRECT, false, false);
 * 
 * // 直接计算（有截断、有PBC）：
 * computeSystemEnergy(state, EnergyMethod::DIRECT, true, true);
 * 
 * // Ewald求和（自动使用PBC和截断）：
 * computeSystemEnergy(state, EnergyMethod::EWALD);
 * 
 * 注意：使用Ewald方法时需要先初始化Ewald参数：
 * initializeEwaldParameters(cutoff, box);
 */

/**
 * @brief 获取能量计算方法的名称
 * 
 * @param method 能量计算方法枚举
 * @return 方法的字符串名称
 */
inline std::string getEnergyMethodName(EnergyMethod method) {
    switch (method) {
        case EnergyMethod::DIRECT: return "Direct";
        case EnergyMethod::EWALD: return "Ewald";
        default: return "Unknown";
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 