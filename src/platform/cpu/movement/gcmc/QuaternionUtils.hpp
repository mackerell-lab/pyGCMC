#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_QUATERNION_UTILS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_QUATERNION_UTILS_HPP

#include "../../../../model/molecular/MCStructures.hpp"
#include <cmath>
#include <random>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

using Quaternion = model::molecular::Quaternion;
using Vector3 = model::molecular::Vector3;

/**
 * @brief Utility functions for quaternion operations in GCMC
 *
 * This class provides quaternion generation and rotation utilities
 * compatible with gcmc_gpu implementation for proper Monte Carlo sampling
 */
class QuaternionUtils {
public:
    /**
     * @brief Generate a uniformly distributed random quaternion
     *
     * Uses the algorithm from Shoemake (1992) "Uniform Random Rotations"
     * Same as gcmc_gpu's create_random_quarternion()
     *
     * @param rng Random number generator
     * @return Random quaternion uniformly distributed on SO(3)
     */
    template<typename RNG>
    static Quaternion generateRandomQuaternion(RNG& rng) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);

        double u = uniform(rng);
        double v = uniform(rng);
        double w = uniform(rng);

        // Using the Shoemake algorithm for uniform distribution
        double sqrt_1_minus_u = std::sqrt(1.0 - u);
        double sqrt_u = std::sqrt(u);
        double two_pi_v = 2.0 * M_PI * v;
        double two_pi_w = 2.0 * M_PI * w;

        // Note: gcmc_gpu returns {x, y, z, w} but we use {w, x, y, z}
        Quaternion q(
            sqrt_u * std::cos(two_pi_w),           // w component
            sqrt_1_minus_u * std::sin(two_pi_v),   // x component
            sqrt_1_minus_u * std::cos(two_pi_v),   // y component
            sqrt_u * std::sin(two_pi_w)            // z component
        );

        q.normalize();
        return q;
    }

    /**
     * @brief Generate a quaternion for small rotation
     *
     * Creates a rotation quaternion with angle uniformly distributed
     * between 0 and maxAngle around a random axis
     *
     * @param rng Random number generator
     * @param maxAngle Maximum rotation angle in radians
     * @return Rotation quaternion
     */
    template<typename RNG>
    static Quaternion generateSmallRotation(RNG& rng, double maxAngle) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        std::normal_distribution<double> normal(0.0, 1.0);

        // Random axis (normalized)
        double ax = normal(rng);
        double ay = normal(rng);
        double az = normal(rng);
        double norm = std::sqrt(ax*ax + ay*ay + az*az);

        if (norm < 1e-10) {
            // Degenerate case, use z-axis
            ax = 0; ay = 0; az = 1; norm = 1;
        }

        ax /= norm;
        ay /= norm;
        az /= norm;

        // Random angle between 0 and maxAngle
        double angle = uniform(rng) * maxAngle;

        // Create quaternion from axis-angle
        double half_angle = angle * 0.5;
        double s = std::sin(half_angle);
        double c = std::cos(half_angle);

        Quaternion q(c, s * ax, s * ay, s * az);
        q.normalize();

        return q;
    }

    /**
     * @brief Create quaternion from Euler angles
     *
     * Compatible with gcmc_gpu's create_quaternion_from_EulerAngle
     *
     * @param yaw Rotation around Z axis (radians)
     * @param pitch Rotation around Y axis (radians)
     * @param roll Rotation around X axis (radians)
     * @return Quaternion representing the rotation
     */
    static Quaternion fromEulerAngles(double yaw, double pitch, double roll) {
        double cy = std::cos(yaw * 0.5);
        double sy = std::sin(yaw * 0.5);
        double cp = std::cos(pitch * 0.5);
        double sp = std::sin(pitch * 0.5);
        double cr = std::cos(roll * 0.5);
        double sr = std::sin(roll * 0.5);

        Quaternion q(
            cr * cp * cy + sr * sp * sy,  // w
            sr * cp * cy - cr * sp * sy,  // x
            cr * sp * cy + sr * cp * sy,  // y
            cr * cp * sy - sr * sp * cy   // z
        );

        q.normalize();
        return q;
    }

    /**
     * @brief Create quaternion from axis and angle
     *
     * Compatible with gcmc_gpu's create_quarternion_from_axis
     *
     * @param axis Rotation axis (will be normalized)
     * @param angle Rotation angle in radians
     * @return Quaternion representing the rotation
     */
    static Quaternion fromAxisAngle(const Vector3& axis, double angle) {
        Vector3 normalized_axis = axis;
        double norm = normalized_axis.norm();

        if (norm < 1e-10) {
            // No rotation
            return Quaternion(1, 0, 0, 0);
        }

        normalized_axis = normalized_axis * (1.0 / norm);

        double half_angle = angle * 0.5;
        double s = std::sin(half_angle);
        double c = std::cos(half_angle);

        Quaternion q(c,
                     s * normalized_axis.x,
                     s * normalized_axis.y,
                     s * normalized_axis.z);

        q.normalize();
        return q;
    }

    /**
     * @brief Apply quaternion rotation to a vector
     *
     * Uses the formula: v' = q * v * q^*
     * Compatible with gcmc_gpu's quarternion_rotation
     *
     * @param q Quaternion
     * @param v Vector to rotate
     * @return Rotated vector
     */
    static Vector3 rotateVector(const Quaternion& q, const Vector3& v) {
        // Using the optimized formula from gcmc_gpu
        double qw = q.w, qx = q.x, qy = q.y, qz = q.z;
        double vx = v.x, vy = v.y, vz = v.z;

        Vector3 result;
        result.x = 2 * vx * (qw*qw + qx*qx - 0.5) +
                   2 * vy * (qx*qy - qw*qz) +
                   2 * vz * (qw*qy + qx*qz);

        result.y = 2 * vx * (qw*qz + qx*qy) +
                   2 * vy * (qw*qw + qy*qy - 0.5) +
                   2 * vz * (qy*qz - qw*qx);

        result.z = 2 * vx * (qx*qz - qw*qy) +
                   2 * vy * (qw*qx + qy*qz) +
                   2 * vz * (qw*qw + qz*qz - 0.5);

        return result;
    }

    /**
     * @brief Check if quaternion represents a valid rotation
     *
     * @param q Quaternion to check
     * @param tolerance Tolerance for unit norm check
     * @return true if quaternion is normalized, false otherwise
     */
    static bool isValidRotation(const Quaternion& q, double tolerance = 1e-6) {
        double norm = q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z;
        return std::abs(norm - 1.0) < tolerance;
    }
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_QUATERNION_UTILS_HPP
