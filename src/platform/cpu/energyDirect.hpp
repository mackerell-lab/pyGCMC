// src/platform/cpu/energyDirect.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
extern const float COULOMB;
extern const float MIN_SAFE_DISTANCE;
extern const float MAX_SAFE_ENERGY;

/**
 * @brief 使用直接计算方法计算系统的非键能量
 * 
 * @param state 系统状态
 * @param use_cutoff 是否使用距离截断
 * @param use_pbc 是否使用周期性边界条件
 */
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

/**
 * @brief 使用直接计算方法计算运动残基的非键能量
 * 
 * @param state 系统状态
 * @param use_cutoff 是否使用距离截断
 * @param use_pbc 是否使用周期性边界条件
 */
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

/**
 * @brief 仅计算系统的范德华能量（带截断）
 * 
 * 此函数仅计算范德华相互作用，不计算静电相互作用。
 * 主要用于与Ewald求和方法配合使用，其中静电相互作用单独处理。
 * 
 * @param state 系统状态
 * @param use_cutoff 是否使用距离截断（一般为true）
 * @param use_pbc 是否使用周期性边界条件
 */
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

// 以下是为了兼容旧接口而保留的函数
void computeMovementEnergy(model::MCState& state);
void computeMovementEnergyCutoff(model::MCState& state);
void computeSystemEnergy(model::MCState& state);
void computeSystemEnergyCutoff(model::MCState& state);
void computeSystemEnergyPBC(model::MCState& state);
void computeSystemEnergyPBCCutoff(model::MCState& state);
void computeSystemVdwEnergyCutoff(model::MCState& state);

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 