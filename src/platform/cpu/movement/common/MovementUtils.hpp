#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_UTILS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_UTILS_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <limits>
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
    
    // Apply rotation to a vector (using verified formula)
    Vector3 rotate(const Vector3& v) const {
        // Standard quaternion rotation formula
        double qw = w, qx = x, qy = y, qz = z;
        double vx = v.x, vy = v.y, vz = v.z;
        
        // Rotation matrix form
        double qw2 = qw * qw;
        double qx2 = qx * qx;
        double qy2 = qy * qy;
        double qz2 = qz * qz;
        
        double rx = vx * (qw2 + qx2 - qy2 - qz2) + 
                   vy * 2.0 * (qx * qy - qw * qz) + 
                   vz * 2.0 * (qx * qz + qw * qy);
                   
        double ry = vx * 2.0 * (qx * qy + qw * qz) + 
                   vy * (qw2 - qx2 + qy2 - qz2) + 
                   vz * 2.0 * (qy * qz - qw * qx);
                   
        double rz = vx * 2.0 * (qx * qz - qw * qy) + 
                   vy * 2.0 * (qy * qz + qw * qx) + 
                   vz * (qw2 - qx2 - qy2 + qz2);
        
        return Vector3(rx, ry, rz);
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
        if (deltaE <= 0) return 1.0;
        double logProb = -beta * deltaE;
        if (logProb < std::log(params.logSpaceMin)) {
            return params.logSpaceMin;
        }
        if (logProb > std::log(params.logSpaceMax)) {
            return params.logSpaceMax;
        }
        return std::exp(logProb);
    }

    static double calculateInsertionProbability(
        int n,
        double deltaE,
        double beta,
        double chemPotential,
        double cavityBias,
        double volumeNm3,
        bool /*useLogSpace*/ = true) {
        AcceptanceLogTerms terms;
        terms.beta = beta;
        terms.betaMu = beta * chemPotential;
        terms.deltaE = deltaE;
        terms.countBefore = n;
        double clampedCavity = std::max(cavityBias, 1e-30);
        terms.logCavity = safeLog(clampedCavity);
        terms.logVolume = safeLog(volumeNm3);
        AcceptanceLogResult res = computeInsertionAcceptance(terms);
        return res.probability;
    }

    static double calculateDeletionProbability(
        int n,
        double deltaE,
        double beta,
        double chemPotential,
        double volumeNm3,
        bool /*useLogSpace*/ = true) {
        AcceptanceLogTerms terms;
        terms.beta = beta;
        terms.betaMu = beta * chemPotential;
        terms.deltaE = deltaE;
        terms.countBefore = n;
        terms.logVolume = safeLog(volumeNm3);
        terms.logCavity = 0.0;
        AcceptanceLogResult res = computeDeletionAcceptance(terms);
        return res.probability;
    }

    static double calculateDeletionProbabilityWithCavity(
        int n,
        double deltaE,
        double beta,
        double chemPotential,
        double cavityBias,
        double volumeNm3,
        bool /*useLogSpace*/ = true) {
        AcceptanceLogTerms terms;
        terms.beta = beta;
        terms.betaMu = beta * chemPotential;
        terms.deltaE = deltaE;
        terms.countBefore = n;
        double clampedCavity = std::max(cavityBias, 1e-30);
        terms.logCavity = safeLog(clampedCavity);
        terms.logVolume = safeLog(volumeNm3);
        AcceptanceLogResult res = computeDeletionAcceptance(terms);
        return res.probability;
    }

    static double calculateInsertionProbabilityWithLambda(
        int n,
        double deltaE,
        double beta,
        double chemPotential,
        double cavityBias,
        double volumeNm3,
        double thermalLambdaNm,
        bool /*useLogSpace*/ = true) {
        AcceptanceLogTerms terms;
        terms.beta = beta;
        terms.betaMu = beta * chemPotential;
        terms.deltaE = deltaE;
        terms.countBefore = n;
        double clampedCavity = std::max(cavityBias, 1e-30);
        terms.logCavity = safeLog(clampedCavity);
        terms.logVolume = safeLog(volumeNm3);
        double lambda = (thermalLambdaNm > 0.0 ? thermalLambdaNm : 1.0);
        terms.logLambda3 = 3.0 * safeLog(lambda);
        AcceptanceLogResult res = computeInsertionAcceptance(terms);
        return res.probability;
    }

    static double calculateDeletionProbabilityWithCavityAndLambda(
        int n,
        double deltaE,
        double beta,
        double chemPotential,
        double cavityBias,
        double volumeNm3,
        double thermalLambdaNm,
        bool /*useLogSpace*/ = true) {
        AcceptanceLogTerms terms;
        terms.beta = beta;
        terms.betaMu = beta * chemPotential;
        terms.deltaE = deltaE;
        terms.countBefore = n;
        double clampedCavity = std::max(cavityBias, 1e-30);
        terms.logCavity = safeLog(clampedCavity);
        terms.logVolume = safeLog(volumeNm3);
        double lambda = (thermalLambdaNm > 0.0 ? thermalLambdaNm : 1.0);
        terms.logLambda3 = 3.0 * safeLog(lambda);
        AcceptanceLogResult res = computeDeletionAcceptance(terms);
        return res.probability;
    }

    static double calculateInsertionProbabilityCBMC(
        int n,
        double beta,
        double chemPotential,
        double volumeNm3,
        double logWnew,
        int Keff,
        double cavityBias) {
        AcceptanceLogTerms terms;
        terms.beta = beta;
        terms.betaMu = beta * chemPotential;
        terms.countBefore = n;
        double clampedCavity = std::max(cavityBias, 1e-30);
        terms.logCavity = safeLog(clampedCavity);
        terms.logVolume = safeLog(volumeNm3);
        terms.cbmcTrials = std::max(Keff, 1);
        terms.logWForward = logWnew;
        AcceptanceLogResult res = computeInsertionAcceptance(terms);
        return res.probability;
    }

private:
    struct AcceptanceLogTerms {
        double logProposalForward = 0.0;
        double logProposalReverse = 0.0;
        double beta = 0.0;
        double betaMu = 0.0;
        double deltaE = 0.0;
        double logVolume = 0.0;
        double logCavity = 0.0;
        int    countBefore = 0;
        int    cbmcTrials = 1;
        double logWForward = 0.0;
        double logWReverse = 0.0;
        double logLambda3 = 0.0;
    };

    struct AcceptanceLogResult {
        double probability = 0.0;
        double logRatio = -std::numeric_limits<double>::infinity();
    };

    static inline double safeLog(double value) {
        return std::log(std::max(value, 1e-30));
    }

    static inline double safeExp(double logValue) {
        if (logValue >= 0.0) return 1.0;
        double lowerBound = -700.0;
        return std::exp(std::max(logValue, lowerBound));
    }

    static AcceptanceLogResult computeInsertionAcceptance(const AcceptanceLogTerms& t) {
        int nPlusOne = t.countBefore + 1;
        double logNplus1 = safeLog(static_cast<double>(std::max(nPlusOne, 1)));
        double logK = safeLog(static_cast<double>(std::max(t.cbmcTrials, 1)));

        double logRatio =
            (t.logProposalForward - t.logProposalReverse)
            - t.beta * t.deltaE
            + t.betaMu
            + (t.logVolume + t.logCavity) - logNplus1
            + (t.logWForward - logK)
            - t.logLambda3;

        AcceptanceLogResult res;
        res.logRatio = logRatio;
        res.probability = safeExp(logRatio);
        return res;
    }

    static AcceptanceLogResult computeDeletionAcceptance(const AcceptanceLogTerms& t) {
        AcceptanceLogResult res;
        if (t.countBefore <= 0) {
            res.probability = 0.0;
            res.logRatio = -std::numeric_limits<double>::infinity();
            return res;
        }

        double logN = safeLog(static_cast<double>(t.countBefore));
        double logK = safeLog(static_cast<double>(std::max(t.cbmcTrials, 1)));

        double logRatio =
            (t.logProposalForward - t.logProposalReverse)
            + t.beta * t.deltaE
            - t.betaMu
            + logN - (t.logVolume + t.logCavity)
            + (logK - t.logWReverse)
            + t.logLambda3;

        res.logRatio = logRatio;
        res.probability = safeExp(logRatio);
        return res;
    }
};

/**
 * Rotation utilities
 */
class RotationUtils {
private:
    static std::mt19937 gen_;
    static std::uniform_real_distribution<> dis_;

public:
    /**
     * Set seed for deterministic quaternion generation
     */
    static void setSeed(unsigned int seed) {
        gen_.seed(seed);
    }

    /**
     * Generate a uniform random quaternion
     */
    static Quaternion generateRandomQuaternion() {
        double u1 = dis_(gen_);
        double u2 = dis_(gen_);
        double u3 = dis_(gen_);
        
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

inline double safeLogProbability(double value) {
    constexpr double kMinProb = 1e-30;
    return std::log(std::max(value, kMinProb));
}

} // namespace utils
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_UTILS_HPP
