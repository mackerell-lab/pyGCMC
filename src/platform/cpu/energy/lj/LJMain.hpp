#pragma once

/**
 * @brief LJ模块统一入口 - Lennard-Jones Potential Module
 * 
 * 本文件聚合LJ模块的所有功能，外部只需包含此文件。
 * 
 * 功能包含：
 * 1. 基础LJ势能计算 (LJPotential)
 * 2. 开关函数支持 (LJSwitch)
 * 3. 安全检查和数值稳定性处理
 * 
 * 典型用法：
 *   #include "lj/LJMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu::lj;
 *   double energy = calcLJEnergyBasic(r2, sigma, eps);
 *   double energy_switch = calcLJEnergyWithSwitching(r2, sigma, eps, info);
 * 
 * @note 这是LJ模块的唯一对外接口，外部模块不应直接包含子模块头文件
 */

// 聚合LJ模块的所有子功能
#include "LJPotential.hpp"
#include "LJSwitch.hpp"

// 提供便捷的命名空间别名（可选）
namespace pygcmc::platform::cpu {
    namespace lj_api = lj;  // 简短别名：lj_api::calcLJEnergyBasic(...)
} 