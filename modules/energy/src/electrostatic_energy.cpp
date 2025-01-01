#include <pygcmc/energy/electrostatic_energy.hpp>
#include <cmath>

namespace pygcmc {

double ElectrostaticEnergy::calculate(const Project& project) {
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

double ElectrostaticEnergy::calculate_atom_energy(const Project& project, const ProjectAtom& atom) {
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

double ElectrostaticEnergy::calculate_pair_energy(const ProjectAtom& atom1, const ProjectAtom& atom2) const {
    // Get charges
    double q1 = atom1.get_charge();
    double q2 = atom2.get_charge();
    
    // Calculate distance with periodic boundary conditions
    double dx = atom1.get_x() - atom2.get_x();
    double dy = atom1.get_y() - atom2.get_y();
    double dz = atom1.get_z() - atom2.get_z();
    
    // Apply minimum image convention
    apply_minimum_image(dx, dy, dz);
    
    double r = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // Calculate Coulomb potential
    // E = k * (q1 * q2) / r
    return k_coulomb * q1 * q2 / r;
}

} // namespace pygcmc 