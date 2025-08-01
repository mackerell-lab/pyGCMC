#pragma once

#include "PGPCore.hpp"
#include <memory>
#include <mutex>

namespace pygcmc {
namespace platform {
namespace cpu {

// 使用智能指针管理全局PGP参数
extern std::unique_ptr<PGPParams> pgp_params_ptr;
extern std::mutex pgp_global_mutex;

// 获取PGP参数的线程安全函数
inline PGPParams& getPGPParams() {
    if (!pgp_params_ptr) {
        std::lock_guard<std::mutex> lock(pgp_global_mutex);
        if (!pgp_params_ptr) {
            pgp_params_ptr = std::make_unique<PGPParams>();
        }
    }
    return *pgp_params_ptr;
}

// 重置PGP参数
inline void resetPGPParamsPtr() {
    std::lock_guard<std::mutex> lock(pgp_global_mutex);
    pgp_params_ptr.reset();
    // 下次调用getPGPParams()时会自动创建新实例
}

// 不使用宏，避免递归调用问题

} // namespace cpu
} // namespace platform
} // namespace pygcmc