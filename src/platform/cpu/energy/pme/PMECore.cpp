#include "PMECore.hpp"
#include "PMEParams.hpp"
#include "PMEBSpline.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Global PME parameters instance
PMEParams pme_params;

// Implementation moved to PMEParams.cpp and PMEBSpline.cpp
// This file now only contains the global PME parameters instance

// <agent-hook:pme_core_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 