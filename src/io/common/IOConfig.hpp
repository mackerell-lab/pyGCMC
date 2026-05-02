#ifndef IO_CONFIG_HPP
#define IO_CONFIG_HPP

#include <iostream>
#include <string>

namespace pygcmc {
namespace io {

// Global configuration for IO operations
struct IOConfig {
    // Control whether to print error messages to stderr
    // Default is false to reduce noise during testing
    static bool verbose_errors;

    // Helper function to conditionally print errors
    static void printError(const std::string& message) {
        if (verbose_errors) {
            std::cerr << message << std::endl;
        }
    }
};

// Initialize static member
inline bool IOConfig::verbose_errors = false;

} // namespace io
} // namespace pygcmc

#endif // IO_CONFIG_HPP
