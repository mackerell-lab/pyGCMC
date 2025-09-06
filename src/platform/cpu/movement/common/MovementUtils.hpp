#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_UTILS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_UTILS_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstdint>
#include "MovementParams.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// 3D Vector type
struct Vector3 {
    double x, y, z;
    
    Vector3() : x(0), y(0), z(0) {}
    Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }
    
    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }
    
    Vector3 operator*(double scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }
    
    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    
    double dot(const Vector3& other) const {
        return x*other.x + y*other.y + z*other.z;
    }
};

// Quaternion for rotations
struct Quaternion {
    double w, x, y, z;
    
    Quaternion() : w(1), x(0), y(0), z(0) {}
    Quaternion(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}
    
    void normalize() {
        double norm = std::sqrt(w*w + x*x + y*y + z*z);
        w /= norm; x /= norm; y /= norm; z /= norm;
    }
    
    // CRITICAL ADDITION: Quaternion multiplication for proper rotation composition
    Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w
        );
    }
};

namespace utils {

/**
 * Log-space utilities for numerical stability
 */
class LogSpaceCalculator {
public:
    /**
     * Calculate log(sum(exp(values))) in a numerically stable way
     */
    static double logSumExp(const std::vector<double>& logValues) {
        if (logValues.empty()) return -std::numeric_limits<double>::infinity();
        
        // Find maximum value for numerical stability
        double maxVal = *std::max_element(logValues.begin(), logValues.end());
        
        // Handle case where all values are -inf
        if (std::isinf(maxVal)) return maxVal;
        
        // Compute sum in a stable way
        double sum = 0.0;
        for (double val : logValues) {
            sum += std::exp(val - maxVal);
        }
        
        return maxVal + std::log(sum);
    }
    
    /**
     * Calculate acceptance probability in log-space
     */
    static double logSpaceAcceptance(double deltaE, double beta, const MovementParams& params) {
        // Always accept if energy decreases
        if (deltaE <= 0) return 1.0;
        
        // Calculate log probability
        double logProb = -beta * deltaE;
        
        // Apply bounds for numerical stability
        if (logProb < std::log(params.logSpaceMin)) {
            return params.logSpaceMin;
        }
        if (logProb > std::log(params.logSpaceMax)) {
            return params.logSpaceMax;
        }
        
        return std::exp(logProb);
    }
    
    /**
     * Calculate insertion acceptance probability
     */
    static double calculateInsertionProbability(
        int n,                      // Current number of molecules
        double deltaE,              // Energy change
        double beta,                // 1/kT
        double chemPotential,       // Chemical potential
        double cavityBias,          // Cavity bias factor
        double volumeNm3,           // System volume in nm^3
        bool useLogSpace = true) {
        
        // Calculate ideal gas concentration using thermodynamic relationship
        // B factor = β*μ + ln(V) where V is volume in appropriate units
        double B = beta * chemPotential + std::log(volumeNm3);
        
        if (useLogSpace) {
            // Log-space calculation
            // For cavity-biased insertion: A_ins includes cavityBias factor
            // This comes from detailed balance with biased proposal
            double logProb = std::log(cavityBias) - std::log(n + 1) + B - beta * deltaE;
            return std::min(1.0, std::exp(logProb));
        } else {
            // Direct calculation
            double prob = cavityBias / (n + 1) * std::exp(B - beta * deltaE);
            return std::min(1.0, prob);
        }
    }
    
    /**
     * Calculate deletion acceptance probability
     */
    static double calculateDeletionProbability(
        int n,                      // Current number of molecules
        double deltaE,              // Energy change
        double beta,                // 1/kT
        double chemPotential,       // Chemical potential
        double volumeNm3,           // System volume in nm^3
        bool useLogSpace = true) {
        
        // Calculate ideal gas concentration using thermodynamic relationship
        // B factor = β*μ + ln(V) where V is volume in appropriate units
        double B = beta * chemPotential + std::log(volumeNm3);
        
        if (useLogSpace) {
            // Log-space calculation
            double logProb = std::log(static_cast<double>(n)) - B - beta * deltaE;
            return std::min(1.0, std::exp(logProb));
        } else {
            // Direct calculation
            double prob = n * std::exp(-B - beta * deltaE);
            return std::min(1.0, prob);
        }
    }
    
    /**
     * Calculate deletion acceptance probability with cavity bias
     * Symmetric to insertion with cavity bias
     */
    static double calculateDeletionProbabilityWithCavity(
        int n,                      // Current number of molecules
        double deltaE,              // Energy change
        double beta,                // 1/kT
        double chemPotential,       // Chemical potential
        double cavityBias,          // Cavity bias factor (probability of selecting this position)
        double volumeNm3,           // System volume in nm^3
        bool useLogSpace = true) {
        
        if (n == 0) {
            return 0.0;
        }
        
        // B factor = β*μ + ln(V)
        double B = beta * chemPotential + std::log(volumeNm3);
        
        if (useLogSpace) {
            // For cavity-biased deletion: A_del includes 1/cavityBias factor
            // This is the reverse of insertion to maintain detailed balance
            double logProb = std::log(static_cast<double>(n)) - std::log(std::max(cavityBias, 1e-30))
                            - B - beta * deltaE;
            return std::min(1.0, std::exp(logProb));
        } else {
            // Direct calculation
            double prob = n / std::max(cavityBias, 1e-30) * std::exp(-B - beta * deltaE);
            return std::min(1.0, prob);
        }
    }
    
    /**
     * Calculate insertion acceptance probability with thermal de Broglie wavelength
     * Includes the Λ³ term for absolute calibration
     */
    static double calculateInsertionProbabilityWithLambda(
        int n,                      // Current number of molecules
        double deltaE,              // Energy change
        double beta,                // 1/kT
        double chemPotential,       // Chemical potential
        double cavityBias,          // Cavity bias factor
        double volumeNm3,           // System volume in nm^3
        double thermalLambdaNm,     // Thermal de Broglie wavelength in nm
        bool useLogSpace = true) {
        
        // B = β*μ + ln(V) - 3*ln(Λ)
        double lambda = (thermalLambdaNm > 0.0 ? thermalLambdaNm : 1.0);
        double B = beta * chemPotential + std::log(volumeNm3) - 3.0 * std::log(lambda);
        
        if (useLogSpace) {
            // For cavity-biased insertion with Lambda: A_ins includes cavityBias
            double logProb = std::log(cavityBias) - std::log(n + 1.0) + B - beta * deltaE;
            return std::min(1.0, std::exp(logProb));
        } else {
            double prob = cavityBias / (n + 1.0) * std::exp(B - beta * deltaE);
            return std::min(1.0, prob);
        }
    }
    
    /**
     * Calculate deletion acceptance probability with cavity bias and thermal wavelength
     * Symmetric to insertion with both cavity bias and Λ³
     */
    static double calculateDeletionProbabilityWithCavityAndLambda(
        int n,                      // Current number of molecules
        double deltaE,              // Energy change
        double beta,                // 1/kT
        double chemPotential,       // Chemical potential
        double cavityBias,          // Cavity bias factor (probability of selecting this position)
        double volumeNm3,           // System volume in nm^3
        double thermalLambdaNm,     // Thermal de Broglie wavelength in nm
        bool useLogSpace = true) {
        
        if (n == 0) {
            return 0.0;
        }
        
        // B = β*μ + ln(V) - 3*ln(Λ)
        double lambda = (thermalLambdaNm > 0.0 ? thermalLambdaNm : 1.0);
        double B = beta * chemPotential + std::log(volumeNm3) - 3.0 * std::log(lambda);
        
        if (useLogSpace) {
            // For cavity-biased deletion with Lambda: A_del includes 1/cavityBias
            double logProb = std::log(static_cast<double>(n)) - std::log(std::max(cavityBias, 1e-30))
                            - B - beta * deltaE;
            return std::min(1.0, std::exp(logProb));
        } else {
            double prob = n / std::max(cavityBias, 1e-30) * std::exp(-B - beta * deltaE);
            return std::min(1.0, prob);
        }
    }
    
    /**
     * Calculate CBMC insertion acceptance probability
     * Uses Rosenbluth weight instead of direct energy change
     */
    static double calculateInsertionProbabilityCBMC(
        int n,                      // Current number of molecules (before insertion)
        double beta,                // 1/kT
        double chemPotential,       // Chemical potential in kJ/mol
        double volumeNm3,           // System volume in nm^3
        double logWnew,             // log(sum(exp(-beta*deltaE_i)))
        int Keff,                   // Effective number of valid trials
        double cavityBias) {        // Cavity bias factor (f_n)
        
        // Calculate ideal gas concentration using thermodynamic relationship
        // B factor = β*μ + ln(V) where V is volume in appropriate units
        const double B = beta * chemPotential + std::log(volumeNm3);
        
        // CBMC acceptance formula in log-space
        const double logA = std::log(cavityBias) - std::log(n + 1.0) 
                          + B + logWnew - std::log(static_cast<double>(Keff));
        
        // Return probability (not log)
        return std::min(1.0, std::exp(logA));
    }
};

/**
 * Rotation utilities
 */
class RotationUtils {
public:
    /**
     * Generate a uniform random quaternion
     */
    static Quaternion generateRandomQuaternion() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(0.0, 1.0);
        
        double u1 = dis(gen);
        double u2 = dis(gen);
        double u3 = dis(gen);
        
        Quaternion q;
        q.w = std::sqrt(1 - u1) * std::sin(2 * M_PI * u2);
        q.x = std::sqrt(1 - u1) * std::cos(2 * M_PI * u2);
        q.y = std::sqrt(u1) * std::sin(2 * M_PI * u3);
        q.z = std::sqrt(u1) * std::cos(2 * M_PI * u3);
        
        q.normalize();
        return q;
    }
    
    /**
     * Convert quaternion to 3x3 rotation matrix
     */
    static void quaternionToMatrix(const Quaternion& q, double matrix[3][3]) {
        double w = q.w, x = q.x, y = q.y, z = q.z;
        
        matrix[0][0] = 1 - 2*y*y - 2*z*z;
        matrix[0][1] = 2*x*y - 2*w*z;
        matrix[0][2] = 2*x*z + 2*w*y;
        
        matrix[1][0] = 2*x*y + 2*w*z;
        matrix[1][1] = 1 - 2*x*x - 2*z*z;
        matrix[1][2] = 2*y*z - 2*w*x;
        
        matrix[2][0] = 2*x*z - 2*w*y;
        matrix[2][1] = 2*y*z + 2*w*x;
        matrix[2][2] = 1 - 2*x*x - 2*y*y;
    }
    
    /**
     * Apply rotation matrix to a vector
     */
    static Vector3 rotateVector(const Vector3& v, const double matrix[3][3]) {
        return Vector3(
            matrix[0][0]*v.x + matrix[0][1]*v.y + matrix[0][2]*v.z,
            matrix[1][0]*v.x + matrix[1][1]*v.y + matrix[1][2]*v.z,
            matrix[2][0]*v.x + matrix[2][1]*v.y + matrix[2][2]*v.z
        );
    }
    
    /**
     * Create quaternion from axis-angle representation
     */
    static Quaternion quaternionFromAxisAngle(const Vector3& axis, double angle) {
        double halfAngle = angle * 0.5;
        double s = std::sin(halfAngle);
        return Quaternion(
            std::cos(halfAngle),  // w
            axis.x * s,            // x
            axis.y * s,            // y
            axis.z * s             // z
        );
    }
};

/**
 * Random number utilities
 */
class RandomUtils {
private:
    static std::mt19937& getGenerator() {
        static std::mt19937 gen;
        static bool initialized = false;
        if (!initialized) {
            // Default initialization with random device
            std::random_device rd;
            gen.seed(rd());
            initialized = true;
        }
        return gen;
    }
    
public:
    /**
     * Set the seed for reproducible random numbers
     * @param seed The seed value (0 = use time-based seed)
     */
    static void setSeed(uint64_t seed) {
        if (seed == 0) {
            // Use time-based seed
            auto now = std::chrono::steady_clock::now().time_since_epoch().count();
            getGenerator().seed(static_cast<unsigned int>(now));
        } else {
            getGenerator().seed(static_cast<unsigned int>(seed));
        }
    }
    
    static double uniform(double min = 0.0, double max = 1.0) {
        std::uniform_real_distribution<> dis(min, max);
        return dis(getGenerator());
    }
    
    static int uniformInt(int min, int max) {
        std::uniform_int_distribution<> dis(min, max);
        return dis(getGenerator());
    }
    
    static bool metropolisAccept(double probability) {
        return uniform() < probability;
    }
    
    static Vector3 randomVector(double maxMagnitude) {
        return Vector3(
            uniform(-maxMagnitude, maxMagnitude),
            uniform(-maxMagnitude, maxMagnitude),
            uniform(-maxMagnitude, maxMagnitude)
        );
    }
};

/**
 * Periodic boundary condition utilities
 */
class PBCUtils {
public:
    static Vector3 applyPBC(const Vector3& position, const Vector3& box) {
        return Vector3(
            position.x - box.x * std::floor(position.x / box.x),
            position.y - box.y * std::floor(position.y / box.y),
            position.z - box.z * std::floor(position.z / box.z)
        );
    }
    
    static void applyPBC(float& x, float& y, float& z, float boxArray[3]) {
        // Apply PBC in-place
        x -= boxArray[0] * std::floor(x / boxArray[0]);
        y -= boxArray[1] * std::floor(y / boxArray[1]);
        z -= boxArray[2] * std::floor(z / boxArray[2]);
    }
    
    static double minimumImageDistance(const Vector3& pos1, const Vector3& pos2, const Vector3& box) {
        Vector3 diff = pos2 - pos1;
        
        // Apply minimum image convention
        if (diff.x > box.x * 0.5) diff.x -= box.x;
        if (diff.x < -box.x * 0.5) diff.x += box.x;
        if (diff.y > box.y * 0.5) diff.y -= box.y;
        if (diff.y < -box.y * 0.5) diff.y += box.y;
        if (diff.z > box.z * 0.5) diff.z -= box.z;
        if (diff.z < -box.z * 0.5) diff.z += box.z;
        
        return diff.norm();
    }
};

} // namespace utils
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_UTILS_HPP