#pragma once
#ifndef PYGCMC_SYSTEM_COMMON_SYSTEMINTERFACE_HPP
#define PYGCMC_SYSTEM_COMMON_SYSTEMINTERFACE_HPP

#include <memory>

namespace pygcmc {
namespace system {
namespace common {

/**
 * @brief System type enumeration
 */
enum class SystemKind {
    MOLECULAR,      // Molecular system for structure/topology handling
    MONTE_CARLO,    // Monte Carlo system for GCMC simulations
    NVT,           // NVT ensemble system (future extension)
    GPU            // GPU-accelerated system (future extension)
};

/**
 * @brief Log level enumeration
 */
enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
};

/**
 * @brief Abstract base interface for all system types
 */
class ISystem {
public:
    virtual ~ISystem() = default;

    /**
     * @brief Initialize the system
     */
    virtual void initialize() = 0;

    /**
     * @brief Get system type
     */
    virtual SystemKind getSystemKind() const = 0;

    /**
     * @brief Check if system is initialized
     */
    virtual bool isInitialized() const = 0;
};

/**
 * @brief Abstract interface for logging systems
 */
class ILogger {
public:
    virtual ~ILogger() = default;

    /**
     * @brief Set verbose mode
     */
    virtual void setVerbose(bool enable) = 0;

    /**
     * @brief Set log level
     */
    virtual void setLogLevel(LogLevel level) = 0;

    /**
     * @brief Get current log level
     */
    virtual LogLevel getLogLevel() const = 0;

    /**
     * @brief Check if verbose mode is enabled
     */
    virtual bool isVerbose() const = 0;
};

/**
 * @brief Abstract interface for molecular systems
 */
class IMolecularSystem : public ISystem {
public:
    virtual ~IMolecularSystem() = default;

    SystemKind getSystemKind() const override { return SystemKind::MOLECULAR; }
};

/**
 * @brief Abstract interface for Monte Carlo systems
 */
class IMonteCarloSystem : public ISystem {
public:
    virtual ~IMonteCarloSystem() = default;

    SystemKind getSystemKind() const override { return SystemKind::MONTE_CARLO; }
};

} // namespace common
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_COMMON_SYSTEMINTERFACE_HPP
