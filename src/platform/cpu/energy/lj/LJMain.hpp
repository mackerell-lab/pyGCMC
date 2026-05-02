#pragma once

/**
 * @brief LJ Module Unified Entry Point - Lennard-Jones Potential Module
 *
 * This file aggregates all functionality of the LJ module, external code only needs to include this file.
 *
 * Features include:
 * 1. Basic LJ potential calculation (LJPotential)
 * 2. Switching function support (LJSwitch)
 * 3. Safety checks and numerical stability handling
 *
 * Typical usage:
 *   #include "lj/LJMain.hpp"
 *
 *   using namespace pygcmc::platform::cpu::lj;
 *   double energy = calcLJEnergyBasic(r2, sigma, eps);
 *   double energy_switch = calcLJEnergyWithSwitching(r2, sigma, eps, info);
 *
 * @note This is the only external interface of the LJ module, external modules should not directly include sub-module header files
 */

// Aggregate all sub-functions of the LJ module
#include "LJPotential.hpp"
#include "LJSwitch.hpp"

