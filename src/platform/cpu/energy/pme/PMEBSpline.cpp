#include "PMEBSpline.hpp"
#include "PMEGlobal.hpp"
#include "PMECore.hpp"
#include "PMEConfig.hpp"
#include "platform/Platform.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Initialize B-splines for PME
 */
void initializePMEBsplines() {
    platform::log(LogLevel::INFO, "Initializing B-splines with order = ", pme_params.splineOrder,
                 " and mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]");

    // Ensure spline order is at least 2
    if (pme_params.splineOrder < 2) {
        platform::log(LogLevel::WARNING, "B-spline order must be at least 2, setting to 2");
        pme_params.splineOrder = 2;
    }

    // Calculate volume
    double boxVolume = pme_params.box[0] * pme_params.box[1] * pme_params.box[2];
    if (boxVolume < 1e-10) {
        platform::log(LogLevel::WARNING, "Box volume near zero, using unit volume for boxfactor");
        boxVolume = 1.0;
    }

    // Find maximum grid size
    int nmax = 0;
    for (int dim = 0; dim < 3; dim++) {
        nmax = (pme_params.meshSize[dim] > nmax) ? pme_params.meshSize[dim] : nmax;
        pme_params.bsplineModuli[dim].resize(pme_params.meshSize[dim]);
    }

    // Initialize B-spline calculation arrays
    std::vector<double> data(pme_params.splineOrder, 0.0);
    std::vector<double> bsplines_data(nmax, 0.0);

    // Initialize with standard B-spline coefficients
    data[0] = 1.0;

    // Calculate B-spline coefficients
    for (int k = 3; k < pme_params.splineOrder; k++) {
        double div = 1.0/(k-1.0);
        data[k-1] = 0.0;
        for (int l = 1; l < (k-1); l++) {
            data[k-l-1] = div*(l*data[k-l-2] + (k-l)*data[k-l-1]);
        }
        data[0] = div*data[0];
    }

    // Calculate final coefficients
    double div = 1.0/(pme_params.splineOrder-1.0);
    data[pme_params.splineOrder-1] = 0.0;
    for (int l = 1; l < (pme_params.splineOrder-1); l++) {
        data[pme_params.splineOrder-l-1] = div*(l*data[pme_params.splineOrder-l-2] + (pme_params.splineOrder-l)*data[pme_params.splineOrder-l-1]);
    }
    data[0] = div*data[0];

    // Initialize bsplines_data
    for (int i = 0; i < nmax; i++) {
        bsplines_data[i] = 0.0;
    }
    for (int i = 1; i <= pme_params.splineOrder; i++) {
        bsplines_data[i] = data[i-1];
    }

    // Calculate B-spline moduli for each dimension
    for (int dim = 0; dim < 3; dim++) {
        int ndata = pme_params.meshSize[dim];
        for (int i = 0; i < ndata; i++) {
            double sc = 0.0, ss = 0.0;
            for (int j = 0; j < ndata; j++) {
                double arg = (2.0*M_PI*i*j)/ndata;
                sc += bsplines_data[j]*cos(arg);
                ss += bsplines_data[j]*sin(arg);
            }
            pme_params.bsplineModuli[dim][i] = sc*sc + ss*ss;
        }

        // Improve numerical stability
        for (int i = 0; i < ndata; i++) {
            if (pme_params.bsplineModuli[dim][i] < 1.0e-7) {
                pme_params.bsplineModuli[dim][i] = (pme_params.bsplineModuli[dim][(i-1+ndata)%ndata] +
                                      pme_params.bsplineModuli[dim][(i+1)%ndata])/2.0;
            }
        }
    }

    // Allocate PME grid
    int totalGridPoints = pme_params.meshSize[0] * pme_params.meshSize[1] * pme_params.meshSize[2];
    pme_params.pmeGrid.resize(totalGridPoints);
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));

    platform::log(LogLevel::INFO, "PME B-splines initialized");
}

/**
 * @brief Set PME parameters
 */
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder, double tolerance) {
    pme_params.alpha = alpha;

    // Ensure grid size is a power of 2
    for (int i = 0; i < 3; i++) {
        platform::log(LogLevel::DEBUG, "Requested mesh dimension ", i, " = ", meshSize[i]);
        if ((meshSize[i] & (meshSize[i] - 1)) != 0) {
            int power = 1;
            while (power < meshSize[i]) {
                power *= 2;
            }
            pme_params.meshSize[i] = power;
            platform::log(LogLevel::WARNING, "PME mesh size must be a power of 2. Adjusting dimension ",
                         i, " from ", meshSize[i], " to ", pme_params.meshSize[i]);
        } else {
            pme_params.meshSize[i] = meshSize[i];
        }
        platform::log(LogLevel::DEBUG, "Effective mesh dimension ", i, " = ", pme_params.meshSize[i]);
    }

    // Set other parameters
    if (splineOrder < 3) {
        platform::log(LogLevel::WARNING, "B-spline order less than 3 may lead to poor accuracy. Setting to 3.");
        pme_params.splineOrder = 3;
    } else if (splineOrder > 6) {
        platform::log(LogLevel::WARNING, "B-spline orders > 6 may be computationally expensive.");
        pme_params.splineOrder = std::min(splineOrder, 10);
    } else {
        pme_params.splineOrder = splineOrder;
    }

    pme_params.tolerance = tolerance;
    pme_params.epsilon_r = 1.0;

    platform::log(LogLevel::INFO, "PME parameters set: alpha = ", alpha,
                 ", mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", pme_params.splineOrder,
                 ", tolerance = ", tolerance);

    // Initialize B-splines
    initializePMEBsplines();
}

/**
 * @brief Auto-adjust PME parameters
 */
void autoAdjustPMEParameters(double error_tolerance, double cutoff_distance, const double box[3]) {
    platform::log(LogLevel::INFO, "Auto-adjusting PME parameters for error tolerance ",
                 error_tolerance, " and cutoff ", cutoff_distance);

    // Determine alpha by targeting the real-space error
    double realSpaceError = error_tolerance / 2.0;
    double alpha = 0.0;

    if (realSpaceError < 1e-6) {
        alpha = 3.5 / cutoff_distance;
    } else {
        double x = -std::log(realSpaceError);
        alpha = std::sqrt(x) / cutoff_distance;

        // Refine with Newton iterations
        for (int i = 0; i < 3; i++) {
            double currentError = std::erfc(alpha * cutoff_distance);
            double derivative = -2.0 / std::sqrt(M_PI) *
                               std::exp(-alpha*alpha*cutoff_distance*cutoff_distance) * cutoff_distance;
            alpha -= (currentError - realSpaceError) / derivative;
        }
    }

    // Limit alpha range
    double minAlpha = 1.0 / cutoff_distance;
    double maxAlpha = 3.0 / cutoff_distance;
    alpha = std::max(minAlpha, std::min(maxAlpha, alpha));

    // Determine mesh dimensions
    int meshSize[3];
    double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
    double reciprocalSpaceError = error_tolerance / 2.0;
    double logTerm = -std::log(reciprocalSpaceError);
    double xterm = alpha * minBoxSize / std::sqrt(logTerm);
    double meshSizeFactor = std::pow(M_PI, 1.0/6.0) * std::pow(6.0 * logTerm, 1.0/3.0) / xterm;

    for (int i = 0; i < 3; i++) {
        double meshScale = box[i] / minBoxSize;
        int size = static_cast<int>(std::ceil(meshSizeFactor * meshScale));

        if (size % 2 != 0) size++;
        if (size < 4) size = 4;

        // Ensure power of 2
        int power = 1;
        while (power < size) power *= 2;
        meshSize[i] = power;
    }

    platform::log(LogLevel::INFO, "Auto-selected PME parameters: alpha = ", alpha,
                 ", mesh size = [", meshSize[0], ",", meshSize[1], ",", meshSize[2], "]");

    // Set parameters
    setPMEParameters(alpha, meshSize, DEFAULT_SPLINE_ORDER, error_tolerance);
    initializePMETables(cutoff_distance);
    pme_params.initialized = true;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
