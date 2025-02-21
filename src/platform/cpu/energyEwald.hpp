// src/platform/cpu/energyEwald.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include "energyDirect.hpp"  // For direct calculation methods

namespace pygcmc {
namespace platform {
namespace cpu {

// Ewald parameters
struct EwaldParams {
    float alpha{1.0f};     // Ewald separation parameter (nm^-1)
    int kmax[3]{6,6,6};    // Maximum reciprocal space wave vectors
    float tolerance{1e-5f}; // Precision control
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