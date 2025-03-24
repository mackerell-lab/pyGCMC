#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include <array>
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Define complex number type alias for consistency
using cmplx = std::complex<double>;

// Add declarations for the CustomFFT namespace functions
namespace CustomFFT {
    void fft3D_forward(cmplx* data, int nx, int ny, int nz);
    void fft3D_backward(cmplx* data, int nx, int ny, int nz);
}

// Constants for PME calculation already defined in energyCommon.hpp
// static const int NUM_TABLE_POINTS = 20000;
// static const double TWO_OVER_SQRT_PI = 2.0/std::sqrt(M_PI);
static const int DEFAULT_SPLINE_ORDER = 4;  // B-spline order for PME 

// PME parameters structure
struct PMEParams {
    double alpha{1.0};     // Ewald separation parameter
    int meshSize[3]{64,64,64}; // Default mesh size for PME
    double tolerance{1e-5f};
    bool initialized{false};
    double cutoff{0.0};    
    double epsilon_r{1.0}; // Relative dielectric constant
    double box[3]{1.0, 1.0, 1.0}; // Box dimensions, default to unit box
    
    // Lookup tables for optimization
    std::vector<double> erfcTable;      
    std::vector<double> ewaldScaleTable;
    double ewaldDX;                     
    double ewaldDXInv;
    double erfcDXInv;

    // B-spline parameters
    int splineOrder{DEFAULT_SPLINE_ORDER};   // B-spline order (typically 4-6)
    std::vector<double> bsplineModuli[3];    // Bspline moduli for the 3 dimensions

    // Reciprocal space grid
    std::vector<std::complex<double>> pmeGrid;
    std::vector<double> pmeCharge;
    
    // Table management methods
    void initializeTables(double cutoff);
    void initializeBsplines();
    void setBox(const double newBox[3]);
    double erfcApprox(double r) const;
    double ewaldScaleApprox(double r) const;
    
    // Error estimation methods
    double estimateRealSpaceError() const {
        return std::erfc(alpha * cutoff);
    }
    
    double estimateReciprocalSpaceError(const double box[3]) const {
        double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
        double minMeshSize = std::min(meshSize[0], std::min(meshSize[1], meshSize[2]));
        double error = (minMeshSize / 2.0) * std::sqrt(alpha * minBoxSize) / 10.0;
        error *= std::exp(-(M_PI * minMeshSize / (2.0 * alpha * minBoxSize)) * 
                         (M_PI * minMeshSize / (2.0 * alpha * minBoxSize)));
        return error;
    }
    
    double estimateTotalError(const double box[3]) const {
        return estimateRealSpaceError() + estimateReciprocalSpaceError(box);
    }
};

// Global PME parameters
extern PMEParams pme_params;

// Function declarations
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder = DEFAULT_SPLINE_ORDER, double tolerance = 1e-5);
void autoAdjustPMEParameters(double error_tolerance, double cutoff_distance, const double box[3]);

/**
 * @brief Initialize PME parameters
 * 
 * @param cutoff Cutoff distance
 * @param box Box dimensions
 * @param alpha Ewald separation parameter (if <=0, automatically calculated)
 * @param meshSize Mesh dimensions for PME (if NULL, automatically calculated)
 * @param splineOrder B-spline order (typically 4-6)
 * @param tolerance Precision control parameter
 */
inline void initializePMEParameters(double cutoff, const double box[3], 
                                  double alpha = 0.0, 
                                  const int* meshSize = nullptr,
                                  int splineOrder = DEFAULT_SPLINE_ORDER,
                                  double tolerance = 1e-5) {
    // Set box dimensions - ensure B-spline initialization uses correct volume
    pme_params.setBox(box);
    
    // If alpha is not specified, calculate the optimal value
    if (alpha <= 0.0) {
        // Auto-adjust parameters includes calling initializeTables and initializeBsplines
        autoAdjustPMEParameters(tolerance, cutoff, box);
    } else {
        // Use the specified parameters
        int mSize[3] = {64, 64, 64}; // Default value
        
        // If mesh size is provided, use it
        if (meshSize != nullptr) {
            mSize[0] = meshSize[0];
            mSize[1] = meshSize[1];
            mSize[2] = meshSize[2];
        }
        
        // Set parameters and initialize tables
        setPMEParameters(alpha, mSize, splineOrder, tolerance);
        pme_params.initializeTables(cutoff);
        pme_params.initializeBsplines();
    }
    
    // Log the final parameters for debugging
    platform::log(LogLevel::INFO, "PME parameters initialized: alpha = ", pme_params.alpha,
                 ", mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", pme_params.splineOrder,
                 ", initialized = ", pme_params.initialized);
}

/**
 * @brief Use PME method to calculate system energy
 * 
 * @param state MC state
 */
void computeSystemEnergyPME(model::MCState& state);

/**
 * @brief Use PME method to calculate energy of moving residue
 * 
 * @param state MC state
 */
void computeMovementEnergyPME(model::MCState& state);

// Internal calculation function declarations
std::pair<double, double> calcPairEnergyPME(double r2, double sigma, double eps, double q1, double q2, const model::MCInfo& info, bool is_excluded = false);
double computeReciprocalPME(model::MCState& state);
double computeSelfEnergyPME(model::MCState& state, bool movement_only);
void computeRealSpacePME(model::MCState& state, bool movement_only, bool store_in_residues = true);

// Helper functions for PME
void spreadChargesOntoGrid(model::MCState& state);
void performFFTForward();
void performFFTBackward();
void computeEnergyFromGrid(double& energy, const double box[3]);
void computeBSplineCoefficients(double fractional, int order, std::vector<double>& coefficients);

// Enable for debug output
extern void setEnergyDebugOutput(bool enable);

} // namespace cpu
} // namespace platform
} // namespace pygcmc
