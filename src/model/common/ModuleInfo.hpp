#pragma once

#ifndef PYGCMC_MODEL_COMMON_MODULE_INFO_HPP
#define PYGCMC_MODEL_COMMON_MODULE_INFO_HPP

#include <string>
#include <sstream>

namespace pygcmc {
namespace model {

/**
 * @brief Module Information and Version
 */
namespace info {
    constexpr const char* VERSION = "2.0.0";
    constexpr const char* BUILD_DATE = __DATE__;
    constexpr const char* DESCRIPTION = "Molecular modeling data structures for GCMC simulation";
    
    constexpr int TOTAL_COMPONENTS = 8;
    constexpr const char* COMPONENTS[] = {
        "common", "atom", "residue", "molecule", 
        "topology", "montecarlo", "param", "structure"
    };
}

/**
 * @brief Get module version and build information
 */
inline std::string get_module_info() {
    std::stringstream ss;
    ss << "Model Module v" << info::VERSION << "\n";
    ss << "Build Date: " << info::BUILD_DATE << "\n";
    ss << "Components: " << info::TOTAL_COMPONENTS << "\n";
    ss << "Description: " << info::DESCRIPTION << "\n";
    return ss.str();
}

} // namespace model
} // namespace pygcmc

/**
 * @brief Module Version Macros
 */
#define PYGCMC_MODEL_VERSION_MAJOR 2
#define PYGCMC_MODEL_VERSION_MINOR 0
#define PYGCMC_MODEL_VERSION_PATCH 0
#define PYGCMC_MODEL_VERSION_STRING "2.0.0"

#endif // PYGCMC_MODEL_COMMON_MODULE_INFO_HPP 