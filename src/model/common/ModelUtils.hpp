#pragma once

#ifndef PYGCMC_MODEL_COMMON_UTILS_HPP
#define PYGCMC_MODEL_COMMON_UTILS_HPP

#include <string>
#include <functional>
#include <cmath>

namespace pygcmc {
namespace model {
namespace common {
namespace utils {

/**
 * @brief Hash combination utility (based on boost::hash_combine)
 */
template<typename T>
inline void hash_combine(std::size_t& seed, const T& v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * @brief Generate hash for a pair of values
 */
template<typename T1, typename T2>
inline std::size_t hash_pair(const T1& v1, const T2& v2) {
    std::size_t seed = 0;
    hash_combine(seed, v1);
    hash_combine(seed, v2);
    return seed;
}

/**
 * @brief Compare floating point numbers with tolerance
 */
inline bool double_equals(double a, double b, double tolerance = 1.0e-9) {
    return std::abs(a - b) < tolerance;
}

/**
 * @brief Check if a double value is finite and valid
 */
inline bool is_valid_double(double value) {
    return std::isfinite(value);
}

/**
 * @brief Safe string comparison (handles null/empty strings)
 */
inline bool string_equals(const std::string& a, const std::string& b) {
    return a == b;
}

/**
 * @brief Case-insensitive string comparison
 */
inline bool string_equals_ci(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    return std::equal(a.begin(), a.end(), b.begin(),
                     [](char a, char b) {
                         return std::tolower(a) == std::tolower(b);
                     });
}

/**
 * @brief Generate unique ID from string components
 */
inline std::string generate_id(const std::string& type, int number, 
                              const std::string& suffix = "") {
    std::string id = type + "_" + std::to_string(number);
    if (!suffix.empty()) {
        id += "_" + suffix;
    }
    return id;
}

/**
 * @brief Trim whitespace from string
 */
inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

/**
 * @brief Convert string to uppercase
 */
inline std::string to_upper(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

/**
 * @brief Convert string to lowercase
 */
inline std::string to_lower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

/**
 * @brief Safe numeric conversion with validation
 */
template<typename T>
inline bool safe_convert(const std::string& str, T& result) {
    try {
        if constexpr (std::is_same_v<T, int>) {
            result = std::stoi(str);
        } else if constexpr (std::is_same_v<T, double>) {
            result = std::stod(str);
        } else if constexpr (std::is_same_v<T, float>) {
            result = std::stof(str);
        }
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace utils
} // namespace common
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_COMMON_UTILS_HPP 