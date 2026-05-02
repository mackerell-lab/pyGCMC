#pragma once

/**
 * @file GCMCSimulationImpl.hpp
 * @brief Implementation headers for GCMC simulations
 */

// Include the actual implementation files
#include "GCMCSimulation.hpp"
#include "GCMCSimulationModular.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace impl {

// Re-export classes for internal use
using GCMCSimulation = ::pygcmc::platform::cpu::simulation::GCMCSimulation;
using GCMCSimulationModular = ::pygcmc::platform::cpu::simulation::GCMCSimulationModular;

} // namespace impl
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc
