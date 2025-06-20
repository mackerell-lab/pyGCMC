#include "LJPotential.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace lj {

double calcLJEnergyBasic(double r2, double sigma, double eps) {
    return calculateLJEnergyNoSwitch<double>(
        r2, sigma, eps, 
        double(LJ_MIN_SAFE_DISTANCE), 
        double(LJ_MAX_SAFE_ENERGY)
    );
}

} // namespace lj
} // namespace cpu
} // namespace platform
} // namespace pygcmc 