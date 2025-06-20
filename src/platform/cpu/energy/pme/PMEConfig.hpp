#ifndef PMEPARAMS_HPP
#define PMEPARAMS_HPP

#include "PMECore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Initialize lookup tables for erfc and scaling functions
 * 
 * @param cutoff Cutoff distance
 */
void initializePMETables(double cutoff);

/**
 * @brief Set box dimensions for PME calculations
 * 
 * @param newBox Box dimensions [3]
 */
void setPMEBox(const double newBox[3]);

/**
 * @brief Approximate erfc function using lookup table
 * 
 * @param r Distance
 * @return Approximate erfc value
 */
double erfcApproximate(double r);

/**
 * @brief Approximate electrostatic scaling function using lookup table
 * 
 * @param r Distance  
 * @return Scaling factor
 */
double ewaldScaleApproximate(double r);

/**
 * @brief Estimate real space error
 * 
 * @return Real space error estimate
 */
double estimatePMERealSpaceError();

/**
 * @brief Estimate reciprocal space error
 * 
 * @param box Box dimensions
 * @return Reciprocal space error estimate
 */
double estimatePMEReciprocalSpaceError(const double box[3]);

/**
 * @brief Estimate total PME error
 * 
 * @param box Box dimensions
 * @return Total error estimate
 */
double estimatePMETotalError(const double box[3]);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PMEPARAMS_HPP 