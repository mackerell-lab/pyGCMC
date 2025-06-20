#pragma once

/**
 * @brief PME模块统一入口 - Particle Mesh Ewald Module
 * 
 * 本文件聚合PME模块的所有功能，外部只需包含此文件。
 * 
 * 功能组件：
 * - PMECore: 核心参数结构和基础函数
 * - PMEConfig: 参数设置和配置管理
 * - PMESetup: 初始化函数
 * - PMEInterface: 高级能量计算接口
 * - PMEFFT: 自定义FFT实现
 * - PMESpline: B样条插值函数
 * - PMEGrid: 网格操作和电荷分布
 * - PMEReal: 实空间能量计算  
 * - PMERecip: 倒空间能量计算
 * - PMESelf: 自能修正
 * - PMEComposite: 统一接口和便捷函数
 * 
 * 典型用法：
 *   #include "pme/PMEMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu;
 *   computeSystemEnergyPME(state);
 *   computeMovementEnergyPME(state);
 * 
 * @note AI Agents功能定位指南:
 * - 能量计算: PMEInterface.hpp -> computeSystemEnergyPME, computeMovementEnergyPME
 * - 组件计算: PMEInterface.hpp -> computeReciprocalPME, computeSelfEnergyPME, computeRealSpacePME
 * - 参数设置: PMESetup.hpp -> setPMEParameters, autoAdjustPMEParameters
 * - 初始化: PMESetup.hpp -> initializePMEParameters, initializePMETables, initializePMEBsplines
 * - 实空间: PMEReal.hpp -> computeRealSpaceEnergy, calcPairEnergyPME
 * - 自能: PMESelf.hpp -> computeSelfEnergyPME, calculateParticleSelfEnergy
 * - 倒空间: PMERecip.hpp -> 倒空间计算
 * - 网格操作: PMEGrid.hpp, PMEGridMap.hpp -> 网格管理和电荷分布
 * - 参数结构: PMECore.hpp -> PMEParams结构和基础函数
 */

// 聚合PME模块的所有子功能
#include "PMECore.hpp"
#include "PMEConfig.hpp"
#include "PMESetup.hpp"
#include "PMEInterface.hpp"
#include "PMEFFT.hpp"
#include "PMESpline.hpp"
#include "PMEGrid.hpp"
#include "PMEReal.hpp"
#include "PMERecip.hpp"
#include "PMESelf.hpp"
#include "PMEComposite.hpp"

// <agent-hook:pme_main>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief PME Module Information
 */
namespace PMEInfo {
    constexpr const char* VERSION = "1.0.0";
    constexpr const char* DESCRIPTION = "Modular Particle Mesh Ewald Implementation";
    constexpr int NUM_MODULES = 9;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "PMECore",
        "PMESetup",
        "PMEInterface",
        "PMEFFT", 
        "PMESpline",
        "PMEGrid",
        "PMEReciprocal",
        "PMERealSpace",
        "PMESelf"
    };
}

} // namespace cpu
} // namespace platform  
} // namespace pygcmc 