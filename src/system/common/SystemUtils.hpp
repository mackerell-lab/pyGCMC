#pragma once
#ifndef PYGCMC_SYSTEM_COMMON_SYSTEMUTILS_HPP
#define PYGCMC_SYSTEM_COMMON_SYSTEMUTILS_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace pygcmc {
namespace system {
namespace common {

/**
 * @brief Utility functions for system operations
 */
class SystemUtils {
public:
    // Periodic boundary condition utilities
    static void applyPBC(float& x, float& y, float& z, const float box[3]);
    static float getMinImageDistSqr(float dx, float dy, float dz, const float box[3]);
    static void wrapCoordinates(float& x, float& y, float& z, const float box[3]);

    // String utilities
    static std::string trim(const std::string& str);
    static std::string toUpper(const std::string& str);
    static std::string toLower(const std::string& str);
    static std::vector<std::string> split(const std::string& str, char delimiter);
    static bool startsWith(const std::string& str, const std::string& prefix);
    static bool endsWith(const std::string& str, const std::string& suffix);

    // Mathematical utilities
    static float distance(float x1, float y1, float z1, float x2, float y2, float z2);
    static float distanceSquared(float x1, float y1, float z1, float x2, float y2, float z2);
    static void normalize(float& x, float& y, float& z);
    static float magnitude(float x, float y, float z);

    // Vector operations
    static float dotProduct(const float a[3], const float b[3]);
    static void crossProduct(const float a[3], const float b[3], float result[3]);

    // Comparison utilities
    static bool isEqual(float a, float b, float tolerance = 1e-6f);
    static bool isZero(float value, float tolerance = 1e-6f);

    // File path utilities
    static std::string getFileExtension(const std::string& filename);
    static std::string getBaseName(const std::string& path);
    static std::string getDirName(const std::string& path);
};

} // namespace common
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_COMMON_SYSTEMUTILS_HPP
