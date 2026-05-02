#include "PMEConfig.hpp"
#include "PMEGlobal.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Initialize lookup tables for erfc and scaling functions
 */
void initializePMETables(double cutoff) {
    pme_params.cutoff = cutoff;
    pme_params.initialized = true;

    // Initialize erfc table for optimization
    pme_params.erfcTable.resize(NUM_TABLE_POINTS);
    pme_params.ewaldScaleTable.resize(NUM_TABLE_POINTS);

    double tableRange = cutoff;
    pme_params.ewaldDX = tableRange / (NUM_TABLE_POINTS - 1);
    pme_params.ewaldDXInv = (NUM_TABLE_POINTS - 1) / tableRange;
    pme_params.erfcDXInv = pme_params.ewaldDXInv;

    // Populate lookup tables
    for (int i = 0; i < NUM_TABLE_POINTS; i++) {
        double r = i * pme_params.ewaldDX;
        double alphaR = pme_params.alpha * r;
        pme_params.erfcTable[i] = std::erfc(alphaR);

        // Store the scale factor for electrostatics
        if (r > 1e-6) {
            pme_params.ewaldScaleTable[i] = pme_params.erfcTable[i] / r;
        } else {
            pme_params.ewaldScaleTable[i] = 2.0 * pme_params.alpha / std::sqrt(M_PI);
        }
    }

    platform::log(LogLevel::INFO, "PME tables initialized with cutoff = ", cutoff,
                 ", alpha = ", pme_params.alpha, ", mesh size = [",
                 pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]");
}

/**
 * @brief Set box dimensions for PME calculations
 */
void setPMEBox(const double newBox[3]) {
    for (int i = 0; i < 3; i++) {
        pme_params.box[i] = newBox[i];
    }
    platform::log(LogLevel::INFO, "PME box dimensions set to [", pme_params.box[0], ", ", pme_params.box[1], ", ", pme_params.box[2], "]");
}

/**
 * @brief Approximate erfc function using lookup table
 */
double erfcApproximate(double r) {
    if (r >= pme_params.cutoff) return 0.0;

    double x = r * pme_params.erfcDXInv;
    int index = static_cast<int>(x);
    if (index >= static_cast<int>(pme_params.erfcTable.size()) - 1) {
        return pme_params.erfcTable.back();
    }

    double fraction = x - index;
    return pme_params.erfcTable[index] + fraction * (pme_params.erfcTable[index+1] - pme_params.erfcTable[index]);
}

/**
 * @brief Approximate electrostatic scaling function using lookup table
 */
double ewaldScaleApproximate(double r) {
    if (r >= pme_params.cutoff) return 0.0;

    double x = r * pme_params.ewaldDXInv;
    int index = static_cast<int>(x);
    if (index >= static_cast<int>(pme_params.ewaldScaleTable.size()) - 1) {
        return pme_params.ewaldScaleTable.back();
    }

    double fraction = x - index;
    return pme_params.ewaldScaleTable[index] + fraction * (pme_params.ewaldScaleTable[index+1] - pme_params.ewaldScaleTable[index]);
}

/**
 * @brief Estimate real space error
 */
double estimatePMERealSpaceError() {
    return std::erfc(pme_params.alpha * pme_params.cutoff);
}

/**
 * @brief Estimate reciprocal space error
 */
double estimatePMEReciprocalSpaceError(const double box[3]) {
    double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
    double minMeshSize = std::min(pme_params.meshSize[0], std::min(pme_params.meshSize[1], pme_params.meshSize[2]));
    double error = (minMeshSize / 2.0) * std::sqrt(pme_params.alpha * minBoxSize) / 10.0;
    error *= std::exp(-(M_PI * minMeshSize / (2.0 * pme_params.alpha * minBoxSize)) *
                     (M_PI * minMeshSize / (2.0 * pme_params.alpha * minBoxSize)));
    return error;
}

/**
 * @brief Estimate total error
 */
double estimatePMETotalError(const double box[3]) {
    return estimatePMERealSpaceError() + estimatePMEReciprocalSpaceError(box);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
