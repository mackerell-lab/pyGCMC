#include "SystemUtils.hpp"
#include <sstream>
#include <cmath>

namespace pygcmc {
namespace system {
namespace common {

// Periodic boundary condition utilities
void SystemUtils::applyPBC(float& x, float& y, float& z, const float box[3]) {
    x = x - box[0] * std::floor(x / box[0]);
    y = y - box[1] * std::floor(y / box[1]);
    z = z - box[2] * std::floor(z / box[2]);
}

float SystemUtils::getMinImageDistSqr(float dx, float dy, float dz, const float box[3]) {
    dx = dx - box[0] * std::round(dx / box[0]);
    dy = dy - box[1] * std::round(dy / box[1]);
    dz = dz - box[2] * std::round(dz / box[2]);
    return dx * dx + dy * dy + dz * dz;
}

void SystemUtils::wrapCoordinates(float& x, float& y, float& z, const float box[3]) {
    while (x < 0.0f) x += box[0];
    while (x >= box[0]) x -= box[0];
    while (y < 0.0f) y += box[1];
    while (y >= box[1]) y -= box[1];
    while (z < 0.0f) z += box[2];
    while (z >= box[2]) z -= box[2];
}

// String utilities
std::string SystemUtils::trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

std::string SystemUtils::toUpper(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

std::string SystemUtils::toLower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

std::vector<std::string> SystemUtils::split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

bool SystemUtils::startsWith(const std::string& str, const std::string& prefix) {
    return str.size() >= prefix.size() &&
           str.compare(0, prefix.size(), prefix) == 0;
}

bool SystemUtils::endsWith(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Mathematical utilities
float SystemUtils::distance(float x1, float y1, float z1, float x2, float y2, float z2) {
    return std::sqrt(distanceSquared(x1, y1, z1, x2, y2, z2));
}

float SystemUtils::distanceSquared(float x1, float y1, float z1, float x2, float y2, float z2) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;
    return dx * dx + dy * dy + dz * dz;
}

void SystemUtils::normalize(float& x, float& y, float& z) {
    float mag = magnitude(x, y, z);
    if (mag > 0.0f) {
        x /= mag;
        y /= mag;
        z /= mag;
    }
}

float SystemUtils::magnitude(float x, float y, float z) {
    return std::sqrt(x * x + y * y + z * z);
}

// Vector operations
float SystemUtils::dotProduct(const float a[3], const float b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void SystemUtils::crossProduct(const float a[3], const float b[3], float result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

// Comparison utilities
bool SystemUtils::isEqual(float a, float b, float tolerance) {
    return std::abs(a - b) < tolerance;
}

bool SystemUtils::isZero(float value, float tolerance) {
    return std::abs(value) < tolerance;
}

// File path utilities
std::string SystemUtils::getFileExtension(const std::string& filename) {
    size_t pos = filename.find_last_of('.');
    if (pos != std::string::npos) {
        return filename.substr(pos + 1);
    }
    return "";
}

std::string SystemUtils::getBaseName(const std::string& path) {
    size_t pos = path.find_last_of("/\\");
    if (pos != std::string::npos) {
        return path.substr(pos + 1);
    }
    return path;
}

std::string SystemUtils::getDirName(const std::string& path) {
    size_t pos = path.find_last_of("/\\");
    if (pos != std::string::npos) {
        return path.substr(0, pos);
    }
    return ".";
}

} // namespace common
} // namespace system
} // namespace pygcmc
