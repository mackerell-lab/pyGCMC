#pragma once

#ifndef PYGCMC_MODEL_COMMON_CONSTANTS_HPP
#define PYGCMC_MODEL_COMMON_CONSTANTS_HPP

#include <limits>

namespace pygcmc {
namespace model {
namespace common {

/**
 * @brief Physical and chemical constants for molecular modeling
 */
namespace constants {

// Physical constants
constexpr double AVOGADRO = 6.02214076e23;           // Avogadro's number (mol^-1)
constexpr double BOLTZMANN = 1.380649e-23;           // Boltzmann constant (J/K)
constexpr double GAS_CONSTANT = 8.314462618;         // Gas constant (J/mol/K)
constexpr double PLANCK = 6.62607015e-34;            // Planck constant (J·s)
constexpr double ELECTRON_CHARGE = 1.602176634e-19;  // Elementary charge (C)

// Unit conversions
constexpr double KCAL_TO_JOULE = 4184.0;             // kcal/mol to J/mol
constexpr double JOULE_TO_KCAL = 1.0 / KCAL_TO_JOULE;
constexpr double ANGSTROM_TO_METER = 1.0e-10;         // Å to m
constexpr double METER_TO_ANGSTROM = 1.0e10;          // m to Å

// Energy unit conversions
constexpr double HARTREE_TO_KCAL = 627.5094740631;    // Hartree to kcal/mol
constexpr double KCAL_TO_HARTREE = 1.0 / HARTREE_TO_KCAL;

// Default values and tolerances
constexpr double DEFAULT_TOLERANCE = 1.0e-6;          // Default numerical tolerance
constexpr double COORDINATE_TOLERANCE = 1.0e-8;       // Coordinate comparison tolerance
constexpr double ENERGY_TOLERANCE = 1.0e-9;           // Energy comparison tolerance

// Invalid/unset values
constexpr double INVALID_DOUBLE = std::numeric_limits<double>::quiet_NaN();
constexpr int INVALID_INT = -1;

// PDB format constants
constexpr double DEFAULT_OCCUPANCY = 1.0;
constexpr double DEFAULT_TEMPFACTOR = 0.0;
constexpr char DEFAULT_ALTLOC = ' ';
constexpr char DEFAULT_CHAIN = ' ';
constexpr char DEFAULT_INSCODE = ' ';

// CHARMM force field constants
constexpr double DEFAULT_WMAIN = 1.0;
constexpr double DEFAULT_WCOMP = 1.0;
constexpr int DEFAULT_MOVE = 1;
constexpr int DEFAULT_IGNORE = 0;
constexpr int DEFAULT_CONSTRAIN = 0;

} // namespace constants
} // namespace common
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_COMMON_CONSTANTS_HPP 