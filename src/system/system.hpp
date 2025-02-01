#pragma once
#ifndef PYGCMC_SYSTEM_SYSTEM_HPP
#define PYGCMC_SYSTEM_SYSTEM_HPP

#include "model/param.hpp"
#include <sstream>
#include <iostream>

namespace pygcmc {
namespace system {

// 日志级别枚举
enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
};

class System {
public:
    // 静态日志控制
    static bool verbose_;
    static LogLevel log_level_;
    
    // 日志功能
    static void set_verbose(bool verbose) { verbose_ = verbose; }
    static void set_log_level(LogLevel level) { log_level_ = level; }
    
    template<typename... Args>
    static void log(LogLevel level, Args... args) {
        if (!verbose_ || level < log_level_) return;
        
        std::stringstream ss;
        (ss << ... << args);
        
        switch (level) {
            case LogLevel::DEBUG:
                std::cout << "[DEBUG] ";
                break;
            case LogLevel::INFO:
                std::cout << "[INFO] ";
                break;
            case LogLevel::WARNING:
                std::cout << "[WARNING] ";
                break;
            case LogLevel::ERROR:
                std::cout << "[ERROR] ";
                break;
        }
        std::cout << ss.str() << std::endl;
    }

    // 从param.hpp移过来的复杂功能
    void initialize_parameters();
    void process_cavity_list();
    void initialize_mc_time_list();

private:
    model::Param params_;  // 系统参数
    // ... 其他现有的成员 ...
};

// 静态成员初始化
inline bool System::verbose_ = false;
inline LogLevel System::log_level_ = LogLevel::INFO;

} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_SYSTEM_HPP
