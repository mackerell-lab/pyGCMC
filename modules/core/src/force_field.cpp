#include "pygcmc/core/force_field.hpp"
#include <cmath>

namespace pygcmc {
namespace core {

double LennardJonesForce::calculate_forces(const System& system,
                                         std::vector<std::vector<std::array<double, 3>>>& forces) const {
    double total_energy = 0.0;

    // Initialize forces array
    for (size_t i = 0; i < system.get_residue_count(); ++i) {
        forces[i].resize(system.get_residue(i).atom_count(), {0.0, 0.0, 0.0});
    }

    // Calculate forces between all pairs of residues
    for (size_t i = 0; i < system.get_residue_count(); ++i) {
        const auto& res1 = system.get_residue(i);
        
        for (size_t j = i + 1; j < system.get_residue_count(); ++j) {
            const auto& res2 = system.get_residue(j);
            
            // Calculate forces between all pairs of particles
            for (size_t a1 = 0; a1 < res1.atom_count(); ++a1) {
                const auto& p1 = res1.particles[a1];
                if (p1.is_virtual) continue;
                
                for (size_t a2 = 0; a2 < res2.atom_count(); ++a2) {
                    const auto& p2 = res2.particles[a2];
                    if (p2.is_virtual) continue;
                    
                    // Compute distance vector
                    auto dr = system.compute_distance(p1, p2);
                    double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    double r = std::sqrt(r2);
                    
                    if (r < 1e-10) continue;  // Skip overlapping particles
                    
                    // Compute LJ force and energy
                    double inv_r = sigma_ / r;
                    double inv_r6 = std::pow(inv_r, 6);
                    double inv_r12 = inv_r6 * inv_r6;
                    
                    double force_scalar = 24.0 * epsilon_ * (2.0 * inv_r12 - inv_r6) / r;
                    double energy = 4.0 * epsilon_ * (inv_r12 - inv_r6);
                    
                    total_energy += energy;
                    
                    // Apply forces
                    for (size_t dim = 0; dim < 3; ++dim) {
                        double force = force_scalar * dr[dim] / r;
                        forces[i][a1][dim] += force;
                        forces[j][a2][dim] -= force;
                    }
                }
            }
        }
    }
    
    return total_energy;
}

double CoulombForce::calculate_forces(const System& system,
                                    std::vector<std::vector<std::array<double, 3>>>& forces) const {
    double total_energy = 0.0;

    // Initialize forces array
    for (size_t i = 0; i < system.get_residue_count(); ++i) {
        forces[i].resize(system.get_residue(i).atom_count(), {0.0, 0.0, 0.0});
    }

    // Calculate forces between all pairs of residues
    for (size_t i = 0; i < system.get_residue_count(); ++i) {
        const auto& res1 = system.get_residue(i);
        
        for (size_t j = i + 1; j < system.get_residue_count(); ++j) {
            const auto& res2 = system.get_residue(j);
            
            // Calculate forces between all pairs of particles
            for (size_t a1 = 0; a1 < res1.atom_count(); ++a1) {
                const auto& p1 = res1.particles[a1];
                if (p1.is_virtual) continue;
                
                for (size_t a2 = 0; a2 < res2.atom_count(); ++a2) {
                    const auto& p2 = res2.particles[a2];
                    if (p2.is_virtual) continue;
                    
                    // Compute distance vector
                    auto dr = system.compute_distance(p1, p2);
                    double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    double r = std::sqrt(r2);
                    
                    if (r < 1e-10) continue;  // Skip overlapping particles
                    
                    // Compute Coulomb force and energy
                    double force_scalar = k_ * p1.charge * p2.charge / (r2 * r);
                    double energy = k_ * p1.charge * p2.charge / r;
                    
                    total_energy += energy;
                    
                    // Apply forces
                    for (size_t dim = 0; dim < 3; ++dim) {
                        double force = force_scalar * dr[dim];
                        forces[i][a1][dim] += force;
                        forces[j][a2][dim] -= force;
                    }
                }
            }
        }
    }
    
    return total_energy;
}

double ResidueBondedForce::calculate_forces(const System& system,
                                          std::vector<std::vector<std::array<double, 3>>>& forces) const {
    double total_energy = 0.0;

    // Initialize forces array
    for (size_t i = 0; i < system.get_residue_count(); ++i) {
        forces[i].resize(system.get_residue(i).atom_count(), {0.0, 0.0, 0.0});
    }

    // Calculate forces for each bond
    for (const auto& bond : bonds_) {
        if (bond.residue1 >= system.get_residue_count() ||
            bond.residue2 >= system.get_residue_count()) {
            continue;
        }

        const auto& res1 = system.get_residue(bond.residue1);
        const auto& res2 = system.get_residue(bond.residue2);

        // Get centers of mass
        auto com1 = res1.center_of_mass();
        auto com2 = res2.center_of_mass();

        // Compute distance vector between centers of mass
        std::array<double, 3> dr;
        for (size_t dim = 0; dim < 3; ++dim) {
            dr[dim] = com2[dim] - com1[dim];
        }

        // Compute distance and force
        double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        double r = std::sqrt(r2);
        
        if (r < 1e-10) continue;  // Skip if residues are too close
        
        double dr_r0 = r - bond.r0;
        double force_scalar = -bond.k * dr_r0;
        double energy = 0.5 * bond.k * dr_r0 * dr_r0;
        
        total_energy += energy;

        // Apply forces to all particles in both residues
        for (size_t dim = 0; dim < 3; ++dim) {
            for (size_t a1 = 0; a1 < res1.atom_count(); ++a1) {
                if (!res1.particles[a1].is_virtual) {
                    forces[bond.residue1][a1][dim] += force_scalar * dr[dim] / r / res1.atom_count();
                }
            }
            
            for (size_t a2 = 0; a2 < res2.atom_count(); ++a2) {
                if (!res2.particles[a2].is_virtual) {
                    forces[bond.residue2][a2][dim] -= force_scalar * dr[dim] / r / res2.atom_count();
                }
            }
        }
    }
    
    return total_energy;
}

double ResidueAngleForce::calculate_forces(const System& system,
                                         std::vector<std::vector<std::array<double, 3>>>& forces) const {
    double total_energy = 0.0;

    // Initialize forces array
    for (size_t i = 0; i < system.get_residue_count(); ++i) {
        forces[i].resize(system.get_residue(i).atom_count(), {0.0, 0.0, 0.0});
    }

    // Calculate forces for each angle
    for (const auto& angle : angles_) {
        if (angle.residue1 >= system.get_residue_count() ||
            angle.residue2 >= system.get_residue_count() ||
            angle.residue3 >= system.get_residue_count()) {
            continue;
        }

        const auto& res1 = system.get_residue(angle.residue1);
        const auto& res2 = system.get_residue(angle.residue2);
        const auto& res3 = system.get_residue(angle.residue3);

        // Get centers of mass
        auto com1 = res1.center_of_mass();
        auto com2 = res2.center_of_mass();
        auto com3 = res3.center_of_mass();

        // Compute vectors between centers of mass
        std::array<double, 3> v1, v2;
        for (size_t dim = 0; dim < 3; ++dim) {
            v1[dim] = com1[dim] - com2[dim];
            v2[dim] = com3[dim] - com2[dim];
        }

        // Compute lengths
        double r1 = std::sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
        double r2 = std::sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
        
        if (r1 < 1e-10 || r2 < 1e-10) continue;  // Skip if residues are too close

        // Compute angle
        double cos_theta = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) / (r1 * r2);
        cos_theta = std::min(1.0, std::max(-1.0, cos_theta));  // Clamp to [-1, 1]
        double theta = std::acos(cos_theta);
        
        // Compute energy
        double dtheta = theta - angle.theta0;
        double energy = 0.5 * angle.k * dtheta * dtheta;
        total_energy += energy;

        // Compute forces
        double sin_theta = std::sin(theta);
        if (std::abs(sin_theta) < 1e-10) continue;  // Skip if angle is 0 or 180 degrees
        
        double force_scalar = -angle.k * dtheta / sin_theta;
        
        // Apply forces to all particles
        for (size_t dim = 0; dim < 3; ++dim) {
            double f1 = force_scalar * (v2[dim]/(r1*r2) - cos_theta*v1[dim]/(r1*r1));
            for (size_t a = 0; a < res1.atom_count(); ++a) {
                if (!res1.particles[a].is_virtual) {
                    forces[angle.residue1][a][dim] += f1 / res1.atom_count();
                }
            }
            
            double f2 = -force_scalar * ((v1[dim]+v2[dim])/(r1*r2) - cos_theta*(v1[dim]/(r1*r1) + v2[dim]/(r2*r2)));
            for (size_t a = 0; a < res2.atom_count(); ++a) {
                if (!res2.particles[a].is_virtual) {
                    forces[angle.residue2][a][dim] += f2 / res2.atom_count();
                }
            }
            
            double f3 = force_scalar * (v1[dim]/(r1*r2) - cos_theta*v2[dim]/(r2*r2));
            for (size_t a = 0; a < res3.atom_count(); ++a) {
                if (!res3.particles[a].is_virtual) {
                    forces[angle.residue3][a][dim] += f3 / res3.atom_count();
                }
            }
        }
    }
    
    return total_energy;
}

} // namespace core
} // namespace pygcmc
