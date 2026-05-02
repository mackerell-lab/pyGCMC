#pragma once

/**
 * @brief Coulomb Module Unified Entry Point - Coulomb Potential Module
 *
 * This file aggregates all functionality of the Coulomb module, external code only needs to include this file.
 *
 * Features include:
 * 1. Basic Coulomb potential calculation (CoulombPotential)
 * 2. Combined LJ and Coulomb pair energy calculation (CoulombPairCore)
 * 3. Safety checks and numerical stability handling
 * 4. Template functions for performance optimization
 *
 * Typical usage:
 *   #include "coulomb/CoulombMain.hpp"
 *
 *   using namespace pygcmc::platform::cpu::coulomb;
 *   double energy = calcCoulombEnergy(r, q1, q2);
 *   auto [vdw, elec] = calcPairEnergy(r2, sigma, eps, q1, q2, info);
 *
 * @note This is the only external interface of the Coulomb module, external modules should not directly include sub-module header files
 */

// Aggregate all sub-functions of the Coulomb module
#include "CoulombPotential.hpp"
#include "CoulombPairCore.hpp"
