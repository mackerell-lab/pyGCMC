#pragma once

#include <memory>
#include <string>
#include <stdexcept>
#include "SystemInterface.hpp"

namespace pygcmc {
namespace system {
namespace common {

/**
 * @brief Factory function for creating system instances
 * 
 * @param kind Type of system to create
 * @return Unique pointer to the created system
 * @throws std::runtime_error if system type is not supported
 */
inline std::unique_ptr<ISystem> createSystem(SystemKind kind) {
    switch (kind) {
        case SystemKind::MOLECULAR:
            // For now, return nullptr as molecular system doesn't implement ISystem yet
            // This can be extended when needed
            throw std::runtime_error("Molecular system factory not yet implemented");
            
        case SystemKind::MONTE_CARLO:
            // For now, return nullptr as Monte Carlo system doesn't implement ISystem yet  
            // This can be extended when needed
            throw std::runtime_error("Monte Carlo system factory not yet implemented");
            
        default:
            throw std::runtime_error("Unknown system kind");
    }
}

/**
 * @brief Get version information for the system module
 * 
 * @return Version string indicating refactored modular design
 */
inline std::string getSystemModuleVersion() {
    return "1.0.0-refactored";
}

} // namespace common
} // namespace system
} // namespace pygcmc 