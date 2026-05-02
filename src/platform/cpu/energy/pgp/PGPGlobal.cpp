#include "PGPGlobal.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Define global variables
std::unique_ptr<PGPParams> pgp_params_ptr;
std::mutex pgp_global_mutex;

} // namespace cpu
} // namespace platform
} // namespace pygcmc
