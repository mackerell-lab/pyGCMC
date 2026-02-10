#include "LogMain.hpp"

namespace pygcmc {
namespace system {
namespace log {

// ILogger interface implementation
void LogMain::setVerbose(bool enable) {
    common::LoggingState::system_enabled = enable;
}

void LogMain::setLogLevel(common::LogLevel level) {
    common::LoggingState::system_level = level;
}

common::LogLevel LogMain::getLogLevel() const {
    return common::LoggingState::system_level;
}

bool LogMain::isVerbose() const {
    return common::LoggingState::system_enabled;
}

// Static interface for backward compatibility
void LogMain::set_verbose(bool verbose) {
    common::LoggingState::system_enabled = verbose;
}

void LogMain::set_log_level(common::LogLevel level) {
    common::LoggingState::system_level = level;
}

bool LogMain::get_verbose() {
    return common::LoggingState::system_enabled;
}

common::LogLevel LogMain::get_log_level() {
    return common::LoggingState::system_level;
}

} // namespace log
} // namespace system
} // namespace pygcmc
