#pragma once

#ifndef PYGCMC_MODEL_CONSTANTS_HPP
#define PYGCMC_MODEL_CONSTANTS_HPP

#include <cmath>
#include <limits>

namespace pygcmc {
namespace model {

/**
 * @brief Physical and chemical constants for GCMC simulations
 */
namespace constants {

    // Physical constants
    constexpr double BOLTZMANN = 0.001987;  ///< Boltzmann constant (kcal/mol/K)
    constexpr double KCAL_TO_KJ = 4.184;    ///< Energy unit conversion (kcal/mol -> kJ/mol)
    constexpr double KJ_TO_KCAL = 1.0 / KCAL_TO_KJ;  ///< Energy unit conversion (kJ/mol -> kcal/mol)
    
    // Conversion factors
    constexpr double ANGSTROM_TO_NM = 0.1;   ///< Length conversion (Å -> nm)
    constexpr double NM_TO_ANGSTROM = 10.0;  ///< Length conversion (nm -> Å)
    constexpr double DEGREE_TO_RADIAN = M_PI / 180.0;  ///< Angle conversion
    constexpr double RADIAN_TO_DEGREE = 180.0 / M_PI;  ///< Angle conversion
    
    // Default values
    constexpr double DEFAULT_TEMPERATURE = 300.0;      ///< Default temperature (K)
    constexpr double DEFAULT_WATER_DENSITY = 55.0;     ///< Water density (M)
    constexpr double DEFAULT_CUTOFF = 12.0;             ///< Default cutoff distance (Å)
    constexpr double DEFAULT_GRID_SPACING = 1.0;       ///< Default grid spacing (Å)
    constexpr double DEFAULT_SIGMA = 2.4;               ///< Default cavity sigma (Å)
    
    // Numerical limits
    constexpr double EPSILON = 1e-8;                    ///< Small value for comparisons
    constexpr double LARGE_NUMBER = 1e30;               ///< Large number
    constexpr double INVALID_VALUE = std::numeric_limits<double>::quiet_NaN();
    
    // CHARMM specific constants
    constexpr int MAX_ATOM_NAME_LENGTH = 4;             ///< Maximum atom name length
    constexpr int MAX_RESIDUE_NAME_LENGTH = 4;          ///< Maximum residue name length
    constexpr int MAX_SEGMENT_NAME_LENGTH = 4;          ///< Maximum segment name length
    
    // MC simulation defaults
    constexpr double DEFAULT_INSERTION_DELETION_FRAC = 0.5;  ///< Default insertion/deletion fraction
    constexpr double DEFAULT_MAX_TRANSLATION = 1.0;          ///< Default max translation (Å)
    constexpr double DEFAULT_MAX_ROTATION = 30.0;            ///< Default max rotation (degrees)
    constexpr unsigned int DEFAULT_CONF_BIAS_TRIALS = 10;    ///< Default configuration bias trials
    
    // Energy calculation defaults
    constexpr double DEFAULT_FRAGMENT_CUTOFF = 10.0;    ///< Default fragment cutoff (Å)
    constexpr double DEFAULT_PROTEIN_CUTOFF = 10.0;     ///< Default protein cutoff (Å)
    constexpr unsigned int DEFAULT_PAIRLIST_FREQ = 1000; ///< Default pairlist update frequency
    
    // File format constants
    constexpr int PDB_ATOM_NAME_WIDTH = 4;               ///< PDB atom name field width
    constexpr int PDB_RESIDUE_NAME_WIDTH = 3;            ///< PDB residue name field width
    constexpr int PDB_CHAIN_WIDTH = 1;                   ///< PDB chain field width
    constexpr int PDB_RESIDUE_NUMBER_WIDTH = 4;          ///< PDB residue number field width
    
    // Validation thresholds
    constexpr double MIN_VALID_MASS = 0.1;               ///< Minimum valid atomic mass
    constexpr double MAX_VALID_MASS = 1000.0;            ///< Maximum valid atomic mass  
    constexpr double MIN_VALID_CHARGE = -10.0;           ///< Minimum valid charge
    constexpr double MAX_VALID_CHARGE = 10.0;            ///< Maximum valid charge
    constexpr double MIN_VALID_COORDINATE = -1e6;        ///< Minimum valid coordinate
    constexpr double MAX_VALID_COORDINATE = 1e6;         ///< Maximum valid coordinate
    
} // namespace constants

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_CONSTANTS_HPP 