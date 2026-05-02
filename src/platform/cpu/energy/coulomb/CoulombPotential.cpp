#include "CoulombPotential.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace coulomb {

double calcCoulombEnergy(double r, double q1, double q2) {
    return calculateCoulombEnergy<double>(r, q1, q2, double(MAX_SAFE_ENERGY));
}

} // namespace coulomb
} // namespace cpu
} // namespace platform
} // namespace pygcmc
