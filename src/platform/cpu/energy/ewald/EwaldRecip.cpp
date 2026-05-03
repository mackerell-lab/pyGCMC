#include "EwaldRecip.hpp"
#include "platform/Platform.hpp"
#include <cmath>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate reciprocal space energy using Ewald summation
 */
double computeReciprocalEnergy(model::MCState& state, bool movement_only) {
    const auto& box = state.info.box;
    const auto& atoms = state.atoms;
    double volume = box[0] * box[1] * box[2];

    // Check system neutrality
    double totalCharge = 0.0;
    for(const auto& atom : atoms) {
        totalCharge += static_cast<double>(atom.charge);
    }
    if (std::abs(totalCharge) > 1e-7) {
        platform::log(LogLevel::WARNING,
            "System charge (", totalCharge, ") is not exactly neutral. ",
            "For better accuracy, consider adjusting charges to ensure strict neutrality.");
    }

    // Constants for reciprocal space calculation
    const double recipCoeff = COULOMB * 4.0 * M_PI / volume;
    const double factorEwald = -1.0 / (4.0 * ewald_params.alpha * ewald_params.alpha);

    double total_energy = 0.0;

    // K-space summation
    for (int rx = -ewald_params.kmax[0]; rx <= ewald_params.kmax[0]; rx++) {
        for (int ry = -ewald_params.kmax[1]; ry <= ewald_params.kmax[1]; ry++) {
            for (int rz = -ewald_params.kmax[2]; rz <= ewald_params.kmax[2]; rz++) {
                // Skip k = 0
                if (rx == 0 && ry == 0 && rz == 0) continue;

                double kx = rx * TWO_PI / box[0];
                double ky = ry * TWO_PI / box[1];
                double kz = rz * TWO_PI / box[2];
                double k2 = kx*kx + ky*ky + kz*kz;

                // Calculate structure factor
                std::complex<double> structureFactor =
                    calculateStructureFactor(state, kx, ky, kz, movement_only);

                double ak = std::exp(k2 * factorEwald) / k2;
                double structureFactorNorm = std::norm(structureFactor);

                // Accumulate energy
                total_energy += recipCoeff * ak * structureFactorNorm;
            }
        }
    }

    // Multiply by 0.5 (standard Ewald convention)
    total_energy *= 0.5;

    return total_energy;
}

/**
 * @brief Calculate structure factor for a given k-vector
 */
std::complex<double> calculateStructureFactor(const model::MCState& state,
                                             double kx, double ky, double kz,
                                             bool movement_only) {
    const auto& atoms = state.atoms;
    int numAtoms = static_cast<int>(atoms.size());

    std::complex<double> structureFactor(0.0, 0.0);

    for (int n = 0; n < numAtoms; n++) {
        if (movement_only) {
            bool in_movement = false;
            for (const auto& movementInfo : state.movementResidues) {
                if (n >= movementInfo.startIndex &&
                    n < movementInfo.startIndex + movementInfo.activeCount) {
                    in_movement = true;
                    break;
                }
            }
            if (!in_movement) continue;
        }

        double kdotr = kx*static_cast<double>(atoms[n].x) +
                       ky*static_cast<double>(atoms[n].y) +
                       kz*static_cast<double>(atoms[n].z);
        std::complex<double> phase(std::cos(kdotr), std::sin(kdotr));
        structureFactor += static_cast<double>(atoms[n].charge) * phase;
    }

    return structureFactor;
}

/**
 * @brief Optimize k-vector summation limits
 */
void optimizeKVectorLimits(const double box[3],
                          double alpha,
                          double tolerance,
                          int kmax[3]) {
    // Calculate optimal kmax for each dimension based on convergence
    double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
    double kmax_float = 2.0 * alpha * minBoxSize *
                       std::sqrt(-std::log(2.0 * tolerance));

    for(int i = 0; i < 3; i++) {
        kmax[i] = static_cast<int>(std::ceil(kmax_float * minBoxSize/box[i]));
        // Ensure minimum value
        if (kmax[i] < 5) kmax[i] = 5;
        // Ensure maximum reasonable value
        if (kmax[i] > 50) kmax[i] = 50;
    }

    platform::log(LogLevel::INFO,
        "Optimized k-vector limits: [", kmax[0], ",", kmax[1], ",", kmax[2], "]");
}

/**
 * @brief Validate reciprocal space calculation parameters
 */
bool validateReciprocalSpaceParameters(const model::MCState& state) {
    // Check if Ewald parameters are initialized
    if (!ewald_params.initialized) {
        platform::log(LogLevel::ERROR, "Ewald parameters not initialized");
        return false;
    }

    // Check box dimensions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        platform::log(LogLevel::ERROR, "Invalid box dimensions for reciprocal space calculation");
        return false;
    }

    // Check k-vector limits
    for(int i = 0; i < 3; i++) {
        if (ewald_params.kmax[i] <= 0) {
            platform::log(LogLevel::ERROR, "Invalid kmax values: must be positive");
            return false;
        }
    }

    // Warn about excessive k-vector count
    long long totalKVectors = (2*ewald_params.kmax[0]+1) *
                             (2*ewald_params.kmax[1]+1) *
                             (2*ewald_params.kmax[2]+1) - 1; // exclude k=0
    if (totalKVectors > 50000) {
        platform::log(LogLevel::WARNING,
            "Large number of k-vectors (", totalKVectors,
            ") may lead to slow reciprocal space calculation");
    }

    return true;
}

/**
 * @brief Get reciprocal space convergence information
 */
void getReciprocalSpaceInfo(const model::MCState& state,
                          int& numKVectors,
                          double& maxKVector,
                          double& estimatedError) {
    const auto& box = state.info.box;

    // Calculate total number of k-vectors
    numKVectors = (2*ewald_params.kmax[0]+1) *
                  (2*ewald_params.kmax[1]+1) *
                  (2*ewald_params.kmax[2]+1) - 1; // exclude k=0

    // Calculate maximum k-vector magnitude
    double kx_max = ewald_params.kmax[0] * TWO_PI / box[0];
    double ky_max = ewald_params.kmax[1] * TWO_PI / box[1];
    double kz_max = ewald_params.kmax[2] * TWO_PI / box[2];
    maxKVector = std::sqrt(kx_max*kx_max + ky_max*ky_max + kz_max*kz_max);

    // Estimate convergence error
    estimatedError = ewald_params.estimateReciprocalSpaceError(
        reinterpret_cast<const double*>(box));

    platform::log(LogLevel::INFO,
        "Reciprocal space info: ", numKVectors, " k-vectors, ",
        "max |k| = ", maxKVector, ", estimated error = ", estimatedError);
}

// <agent-hook:ewald_reciprocal_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
