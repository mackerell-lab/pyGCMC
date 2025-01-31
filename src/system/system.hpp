#pragma once
#ifndef PYGCMC_SYSTEM_SYSTEM_HPP
#define PYGCMC_SYSTEM_SYSTEM_HPP

#include "model/param.hpp"

namespace pygcmc {
namespace system {

class System {
public:
    // ... 其他现有的声明 ...

    // 从param.hpp移过来的复杂功能
    void initialize_parameters();
    void process_cavity_list();
    void initialize_mc_time_list();

private:
    model::Param params_;  // 系统参数
    // ... 其他现有的成员 ...
};

} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_SYSTEM_HPP
