#include "PMEGlobal.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Define global PME parameters smart pointer and mutex
std::unique_ptr<PMEParams> pme_params_ptr;
std::mutex pme_global_mutex;

} // namespace cpu
} // namespace platform
} // namespace pygcmc
