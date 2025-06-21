#include "LogMain.hpp"

namespace pygcmc {
namespace system {
namespace log {

// Static member initialization
bool LogMain::verbose_ = false;
common::LogLevel LogMain::log_level_ = common::LogLevel::INFO;

// ILogger interface implementation
void LogMain::setVerbose(bool enable) {
    instance_verbose_ = enable;
}

void LogMain::setLogLevel(common::LogLevel level) {
    instance_log_level_ = level;
}

common::LogLevel LogMain::getLogLevel() const {
    return instance_log_level_;
}

bool LogMain::isVerbose() const {
    return instance_verbose_;
}

// Static interface for backward compatibility
void LogMain::set_verbose(bool verbose) {
    verbose_ = verbose;
}

void LogMain::set_log_level(common::LogLevel level) {
    log_level_ = level;
}

bool LogMain::get_verbose() {
    return verbose_;
}

common::LogLevel LogMain::get_log_level() {
    return log_level_;
}

} // namespace log
} // namespace system
} // namespace pygcmc 