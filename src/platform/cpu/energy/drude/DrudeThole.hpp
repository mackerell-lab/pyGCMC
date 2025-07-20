#ifndef PYGCMC_PLATFORM_CPU_DRUDE_THOLE_HPP
#define PYGCMC_PLATFORM_CPU_DRUDE_THOLE_HPP

#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * Thole screening functions following OpenMM implementation
 * Reference: OpenMM DrudeForce and Thole (1981) paper
 */
class TholeFunctions {
public:
    /**
     * Calculate Thole screening function and its derivative
     * 
     * @param r Distance between particles
     * @param thole Thole parameter (typically ~1.3)
     * @param alpha1 Polarizability of first particle
     * @param alpha2 Polarizability of second particle
     * @param[out] screening The screening function S(u)
     * @param[out] dScreening_dr The derivative dS/dr
     */
    static void calculateScreening(double r, double thole, double alpha1, double alpha2,
                                   double& screening, double& dScreening_dr) {
        // Calculate screening parameter
        // u = r * a / (alpha1 * alpha2)^(1/6)
        double uscale = thole / std::pow(alpha1 * alpha2, 1.0/6.0);
        double u = r * uscale;
        
        // Linear Thole model: S(u) = 1 - (1 + u/2) * exp(-u)
        double expU = std::exp(-u);
        screening = 1.0 - (1.0 + 0.5 * u) * expU;
        
        // Derivative: dS/du = (1/2) * (1 + u) * exp(-u)
        // dS/dr = (dS/du) * (du/dr) = (dS/du) * uscale
        double dS_du = 0.5 * (1.0 + u) * expU;
        dScreening_dr = dS_du * uscale;
    }
    
    /**
     * Calculate the complete Thole screened interaction energy and forces
     * between two dipoles (4 particle-particle interactions)
     * 
     * This follows the OpenMM convention where:
     * - Drude particle has charge q_drude
     * - Parent atom implicitly has charge -q_drude (offset from its original charge)
     * 
     * The four interactions are:
     * 1. drude1-drude2: +q1*q2 (both negative, product positive)
     * 2. drude1-parent2: -q1*q2 (opposite signs)
     * 3. parent1-drude2: -q1*q2 (opposite signs)
     * 4. parent1-parent2: +q1*q2 (both positive, product positive)
     * 
     * @param drudeCharge1 Charge on first Drude particle
     * @param drudeCharge2 Charge on second Drude particle
     * @param r_dd Distance between Drude particles
     * @param r_dp1 Distance from Drude1 to Parent2
     * @param r_pd1 Distance from Parent1 to Drude2
     * @param r_pp Distance between parent atoms
     * @param thole Thole parameter
     * @param alpha1 Polarizability of first dipole
     * @param alpha2 Polarizability of second dipole
     * @param ONE_4PI_EPS0 Coulomb constant
     * @return Total screened interaction energy
     */
    static double calculateScreenedDipoleInteraction(
        double drudeCharge1, double drudeCharge2,
        double r_dd, double r_dp1, double r_pd1, double r_pp,
        double thole, double alpha1, double alpha2,
        double ONE_4PI_EPS0) {
        
        double energy = 0.0;
        
        // The effective charges for the four interactions
        // Following OpenMM: parent has implicit charge -q_drude
        double q_squared = drudeCharge1 * drudeCharge2;
        
        // Calculate screening for each pair
        double screening, dScreening_dr;
        
        // 1. Drude1-Drude2 (+q1*q2)
        if (r_dd > 1e-10) {
            calculateScreening(r_dd, thole, alpha1, alpha2, screening, dScreening_dr);
            energy += ONE_4PI_EPS0 * q_squared * screening / r_dd;
        }
        
        // 2. Drude1-Parent2 (-q1*q2)
        if (r_dp1 > 1e-10) {
            calculateScreening(r_dp1, thole, alpha1, alpha2, screening, dScreening_dr);
            energy -= ONE_4PI_EPS0 * q_squared * screening / r_dp1;
        }
        
        // 3. Parent1-Drude2 (-q1*q2)
        if (r_pd1 > 1e-10) {
            calculateScreening(r_pd1, thole, alpha1, alpha2, screening, dScreening_dr);
            energy -= ONE_4PI_EPS0 * q_squared * screening / r_pd1;
        }
        
        // 4. Parent1-Parent2 (+q1*q2)
        if (r_pp > 1e-10) {
            calculateScreening(r_pp, thole, alpha1, alpha2, screening, dScreening_dr);
            energy += ONE_4PI_EPS0 * q_squared * screening / r_pp;
        }
        
        return energy;
    }
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_DRUDE_THOLE_HPP