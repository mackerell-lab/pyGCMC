#include "SystemModule.hpp"
#include <stdexcept>

namespace pygcmc {
namespace system {

std::unique_ptr<common::ISystem> createSystem(common::SystemKind kind) {
    switch (kind) {
        case common::SystemKind::MOLECULAR:
            // For now, return nullptr as molecular system doesn't implement ISystem yet
            // This can be extended when needed
            throw std::runtime_error("Molecular system factory not yet implemented");
            
        case common::SystemKind::MONTE_CARLO:
            // For now, return nullptr as Monte Carlo system doesn't implement ISystem yet  
            // This can be extended when needed
            throw std::runtime_error("Monte Carlo system factory not yet implemented");
            
        default:
            throw std::runtime_error("Unknown system kind");
    }
}

void initializeLogging(bool verbose, common::LogLevel level) {
    log::LogMain::set_verbose(verbose);
    log::LogMain::set_log_level(level);
}

std::string getSystemModuleVersion() {
    return "1.0.0-refactored";
}

} // namespace system
} // namespace pygcmc 