#pragma once

/**
 * @brief Ewald模块统一入口 - Ewald Summation Module
 * 
 * 本文件聚合Ewald模块的所有功能，外部只需包含此文件。
 * 
 * 功能组件：
 * - EwaldCore: 核心参数和初始化
 * - EwaldReal: 实空间能量计算
 * - EwaldRecip: 倒空间能量计算  
 * - EwaldSelf: 自能修正
 * - EwaldComposite: 高级统一接口
 * 
 * 典型用法：
 *   #include "ewald/EwaldMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu;
 *   computeSystemEnergyEwald(state);
 *   computeMovementEnergyEwald(state);
 * 
 * @note AI Agents功能定位指南:
 * - 能量计算: EwaldComposite.hpp -> computeSystemEnergyEwald, computeMovementEnergyEwald
 * - 实空间: EwaldReal.hpp -> computeRealSpaceEwald, calcPairEnergyEwaldRealSpace
 * - 倒空间: EwaldRecip.hpp -> computeReciprocalEnergy
 * - 自能: EwaldSelf.hpp -> computeSelfEnergy
 * - 初始化: EwaldInterface.hpp -> initializeEwald, isEwaldInitialized
 * - 参数设置: EwaldCore.hpp -> setEwaldParameters, autoAdjustParameters
 */

// 聚合Ewald模块的所有子功能
#include "EwaldCore.hpp"
#include "EwaldReal.hpp" 
#include "EwaldRecip.hpp"
#include "EwaldSelf.hpp"
#include "EwaldComposite.hpp"

// <agent-hook:ewald_main>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Ewald Module Information
 */
namespace EwaldInfo {
    constexpr const char* VERSION = "1.0.0";
    constexpr const char* DESCRIPTION = "Modular Ewald Summation Implementation";
    constexpr int NUM_MODULES = 4;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "EwaldCore",
        "EwaldRealSpace",
        "EwaldReciprocal", 
        "EwaldSelf"
    };
}

} // namespace cpu
} // namespace platform  
} // namespace pygcmc 