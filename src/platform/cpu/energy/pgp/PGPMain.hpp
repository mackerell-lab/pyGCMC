#pragma once

/**
 * @brief PGP模块统一入口 - Precomputed Grid-Potential Module
 * 
 * 本文件聚合PGP模块的所有功能，外部只需包含此文件。
 * 
 * 功能组件：
 * - PGPCore: 核心参数和初始化
 * - PGPGrid: 网格操作和势能网格
 * - PGPInterpolation: 网格插值和能量计算
 * - PGPPrecompute: 预计算算法优化
 * - PGPReal: 实空间能量计算
 * - PGPSelf: 自能修正
 * - PGPSystem: 完整系统能量评估
 * - PGPComposite: 高级统一接口
 * 
 * 典型用法：
 *   #include "pgp/PGPMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu;
 *   computeSystemEnergyPGP(state);
 *   computeMovementEnergyPGP(state);
 * 
 * @note AI Agents功能定位指南:
 * - 能量计算: PGPSystem.hpp -> computeSystemEnergyPGP, computeMovementEnergyPGP
 * - 网格插值: PGPSystem.hpp -> interpolateMoleculeEnergy, calculateMoleculeEnergy
 * - 网格预计算: PGPPrecompute.hpp -> precomputeGridPotential, setPGPParameters
 * - 实空间: PGPReal.hpp -> computeRealSpacePGP
 * - 自能: PGPSelf.hpp -> computeSelfEnergyPGP
 * - 网格操作: PGPGrid.hpp -> initializePotentialGrid
 * - 参数设置: PGPCore.hpp -> setPGPParameters, PGPParams结构
 */

// 聚合PGP模块的所有子功能
#include "PGPCore.hpp"
#include "PGPGrid.hpp"
#include "PGPInterpolation.hpp"
#include "PGPPrecompute.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"
#include "PGPSystem.hpp"
#include "PGPComposite.hpp"

// <agent-hook:pgp_main>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief PGP Module Information
 */
namespace PGPInfo {
    constexpr const char* VERSION = "1.0.0";
    constexpr const char* DESCRIPTION = "Modular Precomputed Grid-Potential Implementation";
    constexpr int NUM_MODULES = 7;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "PGPCore",
        "PGPGrid",
        "PGPInterpolation", 
        "PGPPrecompute",
        "PGPRealSpace",
        "PGPSelfEnergy",
        "PGPSystemEnergy"
    };
}

} // namespace cpu
} // namespace platform  
} // namespace pygcmc 