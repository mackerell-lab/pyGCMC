#include "PGPGlobal.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// 定义全局变量
std::unique_ptr<PGPParams> pgp_params_ptr;
std::mutex pgp_global_mutex;

} // namespace cpu
} // namespace platform
} // namespace pygcmc