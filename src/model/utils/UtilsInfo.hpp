#pragma once

#ifndef PYGCMC_MODEL_UTILS_INFO_HPP
#define PYGCMC_MODEL_UTILS_INFO_HPP

#include <string>
#include <sstream>
#include <tuple>

namespace pygcmc {
namespace model {

/**
 * @brief Module Information and Version
 */
namespace info {
    constexpr const char* VERSION = "2.0.0";
    constexpr const char* BUILD_DATE = __DATE__;
    constexpr const char* BUILD_TIME = __TIME__;
    constexpr const char* DESCRIPTION = "Molecular modeling data structures for GCMC simulation";
    constexpr const char* MODULE_NAME = "Model Utils";
    
    constexpr int TOTAL_COMPONENTS = 8;
    constexpr const char* COMPONENTS[] = {
        "common", "atom", "residue", "molecule", 
        "topology", "montecarlo", "param", "structure"
    };
}

/**
 * @brief Get model module version
 */
inline std::string getModelVersion() {
    return info::VERSION;
}

/**
 * @brief Get module version and build information
 */
inline std::string get_module_info() {
    std::stringstream ss;
    ss << info::MODULE_NAME << " v" << info::VERSION << "\n";
    ss << "Build Date: " << info::BUILD_DATE << " " << info::BUILD_TIME << "\n";
    ss << "Components: " << info::TOTAL_COMPONENTS << "\n";
    ss << "Description: " << info::DESCRIPTION << "\n";
    return ss.str();
}

/**
 * @brief Get version as tuple (major, minor, patch)
 */
inline std::tuple<int, int, int> get_version_tuple() {
    return std::make_tuple(2, 0, 0);
}

/**
 * @brief Check if version meets requirements
 */
inline bool is_version_at_least(int major_required, int minor_required = 0, int patch_required = 0) {
    auto [major, minor, patch] = get_version_tuple();
    
    if (major > major_required) return true;
    if (major < major_required) return false;
    
    if (minor > minor_required) return true;
    if (minor < minor_required) return false;
    
    return patch >= patch_required;
}

} // namespace model

// Export to parent namespace for convenience
using model::getModelVersion;

} // namespace pygcmc

/**
 * @brief Module Version Macros
 */
#define PYGCMC_MODEL_VERSION_MAJOR 2
#define PYGCMC_MODEL_VERSION_MINOR 0
#define PYGCMC_MODEL_VERSION_PATCH 0
#define PYGCMC_MODEL_VERSION_STRING "2.0.0"

#endif // PYGCMC_MODEL_UTILS_INFO_HPP 