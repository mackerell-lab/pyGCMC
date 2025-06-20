#include "PMESpline.hpp"
#include "PMECore.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute B-spline coefficients exactly like pme.cpp
 * 
 * @param fractional Fractional position (0-1)
 * @param order Spline order
 * @param coefficients Output coefficients
 */
void computeBSplineCoefficients(double fractional, int order, std::vector<double>& coefficients) {
    // Ensure coefficients vector is of correct size
    coefficients.resize(order);
    
    // Zero out all coefficients
    for (int i = 0; i < order; i++) {
        coefficients[i] = 0.0;
    }
    
    // Get fractional part
    double dr = fractional;
    
    // Initialize second-order B-spline basis coefficients
    coefficients[0] = 1.0 - dr;
    coefficients[1] = dr;
    
    // Recursively compute B-spline coefficients from third order to order-1 (excluding last step)
    for (int k = 3; k < order; k++) {
        double div = 1.0 / (k - 1.0);
        coefficients[k-1] = div * dr * coefficients[k-2];
        
        for (int i = 1; i < (k-1); i++) {
            coefficients[k-i-1] = div * ((dr+i) * coefficients[k-i-2] + 
                                         (k-i-dr) * coefficients[k-i-1]);
        }
        
        coefficients[0] = div * (1.0-dr) * coefficients[0];
    }
    
    // Last step: special handling for k=order case
    double div = 1.0 / (order - 1);
    coefficients[order-1] = div * dr * coefficients[order-2];
    
    for (int i = 1; i < (order-1); i++) {
        coefficients[order-i-1] = div * ((dr+i) * coefficients[order-i-2] + 
                                        (order-i-dr) * coefficients[order-i-1]);
    }
    coefficients[0] = div * (1.0-dr) * coefficients[0];
    
    // Verify coefficient sum
    double sum = 0.0;
    for (int i = 0; i < order; i++) {
        sum += coefficients[i];
    }
    
    // Warn only if deviation is significant
    if (std::abs(sum - 1.0) > 1e-5) {
        platform::log(LogLevel::WARNING, "Warning: B-spline coefficient sum (" + std::to_string(sum) + ") deviates significantly from 1");
    }
}

/**
 * @brief Initialize B-spline moduli for PME calculations
 */
void initializeBSplineModuli(PMEParams& params) {
    // This functionality is currently implemented in PMEParams::initializeBsplines()
    // We keep this as a wrapper for consistency
    params.initializeBsplines();
}

/**
 * @brief Compute B-spline interpolation weights for charge spreading
 */
void computeInterpolationWeights(const double position[3], int order, double* weights) {
    // Compute B-spline coefficients for each dimension
    std::vector<double> coefficients(order);
    
    for (int dim = 0; dim < 3; dim++) {
        computeBSplineCoefficients(position[dim], order, coefficients);
        for (int i = 0; i < order; i++) {
            weights[dim * order + i] = coefficients[i];
        }
    }
}

/**
 * @brief Validate B-spline coefficients for consistency
 */
bool validateBSplineCoefficients(const std::vector<double>& coefficients, int order) {
    if (static_cast<int>(coefficients.size()) != order) {
        return false;
    }
    
    double sum = 0.0;
    for (int i = 0; i < order; i++) {
        sum += coefficients[i];
    }
    
    // Check if sum is close to 1.0 (within tolerance)
    return std::abs(sum - 1.0) < 1e-5;
}

// <agent-hook:spline_implementation>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 