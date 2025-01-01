#include <pygcmc/energy/lennard_jones_energy.hpp>
#include <cmath>

namespace pygcmc {

double LennardJonesEnergy::calculate(const Project& project) {
    double total_energy = 0.0;
    const auto& atoms = project.get_atoms();
    
    // Calculate pairwise interactions
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
            total_energy += calculate_pair_energy(atoms[i], atoms[j]);
        }
    }
    
    return total_energy;
}

double LennardJonesEnergy::calculate_atom_energy(const Project& project, const ProjectAtom& atom) {
    double energy = 0.0;
    const auto& atoms = project.get_atoms();
    
    // Calculate interaction with all other atoms
    for (const auto& other_atom : atoms) {
        if (&other_atom != &atom) {
            energy += calculate_pair_energy(atom, other_atom);
        }
    }
    
    return energy;
}

double LennardJonesEnergy::calculate_pair_energy(const ProjectAtom& atom1, const ProjectAtom& atom2) const {
    // Get LJ parameters from force field
    double epsilon1 = atom1.get_epsilon();
    double sigma1 = atom1.get_sigma();
    double epsilon2 = atom2.get_epsilon();
    double sigma2 = atom2.get_sigma();
    
    // Lorentz-Berthelot mixing rules
    double epsilon = std::sqrt(epsilon1 * epsilon2);
    double sigma = (sigma1 + sigma2) * 0.5;
    
    // Calculate distance with periodic boundary conditions
    double dx = atom1.get_x() - atom2.get_x();
    double dy = atom1.get_y() - atom2.get_y();
    double dz = atom1.get_z() - atom2.get_z();
    
    // Apply minimum image convention
    apply_minimum_image(dx, dy, dz);
    
    double r2 = dx*dx + dy*dy + dz*dz;
    
    // Calculate LJ potential
    double sig_r6 = std::pow(sigma*sigma/r2, 3);
    double sig_r12 = sig_r6 * sig_r6;
    
    return 4.0 * epsilon * (sig_r12 - sig_r6);
}

} // namespace pygcmc 