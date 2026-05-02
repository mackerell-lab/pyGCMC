#include "EwaldCore.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

// Global Ewald parameters instance
EwaldParams ewald_params;

// Constants
const double TWO_PI = 2.0 * M_PI;
const double SQRT_PI = std::sqrt(M_PI);

/**
 * @brief Initialize lookup tables for erfc function optimization
 */
void EwaldParams::initializeTables(double cutoff) {
    this->cutoff = cutoff;
    ewaldDX = cutoff/NUM_TABLE_POINTS;
    ewaldDXInv = 1.0/ewaldDX;
    erfcDXInv = 1.0/(ewaldDX*alpha);

    // Debug output
    platform::log(LogLevel::INFO,
        "EwaldParams::initializeTables: cutoff=", cutoff,
        " alpha=", alpha,
        " ewaldDX=", ewaldDX,
        " ewaldDXInv=", ewaldDXInv,
        " erfcDXInv=", erfcDXInv,
        " NUM_TABLE_POINTS=", NUM_TABLE_POINTS);

    erfcTable.resize(NUM_TABLE_POINTS + 4);
    ewaldScaleTable.resize(NUM_TABLE_POINTS + 4);

    // Build erfc lookup table
    for(int i = 0; i < NUM_TABLE_POINTS + 4; i++) {
        double r = i * ewaldDX;
        double alphaR = alpha * r;
        erfcTable[i] = std::erfc(alphaR);

        // Print first 5 and last 5 values for debugging
        if (i < 5 || i > NUM_TABLE_POINTS - 1) {
            platform::log(LogLevel::INFO,
                "erfcTable[", i, "]: r=", r,
                " alphaR=", alphaR,
                " erfc=", erfcTable[i]);
        }
    }
}

/**
 * @brief Initialize exp(ikr) lookup tables for reciprocal space optimization
 */
void EwaldParams::initializeExpIkrTable(int numAtoms) {
    maxK = std::max(kmax[0], std::max(kmax[1], kmax[2]));
    expIkrTable.resize(maxK * numAtoms * 3);
    expIkrXY.resize(numAtoms);
}

/**
 * @brief Fast erfc approximation using direct calculation
 */
double EwaldParams::erfcApprox(double r) const {
    double alphaR = alpha * r;
    double result = std::erfc(alphaR);

    platform::log(LogLevel::DEBUG,
        "EwaldParams::erfcApprox: r=", r,
        " alpha=", alpha,
        " alpha*r=", alphaR,
        " erfc=", result);

    return result;
}

/**
 * @brief Ewald scaling factor approximation using lookup table
 */
double EwaldParams::ewaldScaleApprox(double r) const {
    double x = r * ewaldDXInv;
    int index = std::min(static_cast<int>(x), NUM_TABLE_POINTS);
    double coeff2 = x - index;
    double coeff1 = 1.0 - coeff2;
    return coeff1 * ewaldScaleTable[index] + coeff2 * ewaldScaleTable[index + 1];
}

/**
 * @brief Estimate real space error
 */
double EwaldParams::estimateRealSpaceError() const {
    return std::erfc(alpha * cutoff);
}

/**
 * @brief Estimate reciprocal space error
 */
double EwaldParams::estimateReciprocalSpaceError(const double box[3]) const {
    double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
    double minKmax = std::min(kmax[0], std::min(kmax[1], kmax[2]));
    double error = minKmax * std::sqrt(alpha * minBoxSize) / 20.0;
    error *= std::exp(-(M_PI * minKmax / (alpha * minBoxSize)) *
                     (M_PI * minKmax / (alpha * minBoxSize)));
    return error;
}

/**
 * @brief Estimate total error
 */
double EwaldParams::estimateTotalError(const double box[3]) const {
    return estimateRealSpaceError() + estimateReciprocalSpaceError(box);
}

/**
 * @brief Auto-adjust Ewald parameters for optimal performance
 */
void autoAdjustParameters(double error_tolerance, double cutoff_distance, const double box[3]) {
    // Check if cutoff is less than half the box length
    double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
    if (cutoff_distance >= 0.5 * minBoxSize) {
        throw std::runtime_error("Cutoff distance must be less than half the smallest box dimension");
    }

    // Calculate optimal alpha based on error tolerance and cutoff
    ewald_params.alpha = std::sqrt(-std::log(2.0 * error_tolerance)) / cutoff_distance;

    // Calculate optimal kmax for each dimension
    double kmax_float = 2.0 * ewald_params.alpha * minBoxSize *
                       std::sqrt(-std::log(2.0 * error_tolerance));

    for(int i = 0; i < 3; i++) {
        ewald_params.kmax[i] = static_cast<int>(std::ceil(kmax_float * minBoxSize/box[i]));
    }

    // Initialize lookup tables
    ewald_params.initializeTables(cutoff_distance);
    ewald_params.initialized = true;
}

/**
 * @brief Set Ewald calculation parameters explicitly
 */
void setEwaldParameters(double alpha, const int kmax[3], double tolerance) {
    ewald_params.alpha = alpha;
    for(int i = 0; i < 3; i++) {
        ewald_params.kmax[i] = kmax[i];
    }
    ewald_params.tolerance = tolerance;

    // Re-initialize real-space lookup tables with the new alpha
    if (ewald_params.cutoff <= 0.0) {
        ewald_params.cutoff = 1.2;  // Default 1.2 nm cutoff
    }
    ewald_params.initializeTables(ewald_params.cutoff);

    ewald_params.initialized = true;
}

/**
 * @brief Initialize Ewald parameters wrapper function
 */
void initializeEwaldParameters(double cutoff, const double box[3],
                             double alpha, double tolerance) {
    if (alpha <= 0.0) {
        autoAdjustParameters(tolerance, cutoff, box);
    } else {
        int kmax[3] = {15, 15, 15}; // Default value
        setEwaldParameters(alpha, kmax, tolerance);
        ewald_params.initializeTables(cutoff);
    }
}

// <agent-hook:ewald_core_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
