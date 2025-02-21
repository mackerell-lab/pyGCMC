// src/platform/cpu/energyEwald.hpp
#pragma once

#include "../../model/montecarlo.hpp"
#include "../platform.hpp"
#include "energyCommon.hpp"
#include "energyDirect.hpp"  // For direct calculation methods

namespace pygcmc {
namespace platform {
namespace cpu {

// Ewald parameters
struct EwaldParams {
    float alpha{1.0f};     // Ewald分离参数 (nm^-1)
    int kmax[3]{6,6,6};    // 倒空间最大波矢
    float tolerance{1e-5f}; // 精度控制
    bool initialized{false};
};

// Global Ewald parameters
extern EwaldParams ewald_params;

// Function to set Ewald parameters
void setEwaldParameters(float alpha, const int kmax[3], float tolerance = 1e-5f);

// Ewald method interfaces
void computeSystemEnergyEwald(model::MCState& state);
void computeMovementEnergyEwald(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 