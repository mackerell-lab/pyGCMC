// modules/core/src/structure.cpp

#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"
#include <memory>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace pygcmc {
namespace core {

void Structure::apply_forcefield(const std::shared_ptr<ForceField>& forcefield) {
    if (!forcefield) {
        throw std::runtime_error("Cannot apply null forcefield");
    }
    
    // Store the forcefield
    forcefield_ = forcefield;
    
    // Apply forcefield parameters to each atom
    for (auto& atom : atoms_) {
        if (atom && !atom->topo_type.empty()) {
            // Get parameters for this atom type
            const auto& params = forcefield->nonbonded_params();
            auto it = params.find(atom->topo_type);
            if (it != params.end()) {
                atom->forcefield_epsilon = it->second.epsilon;
                atom->forcefield_rmin = it->second.rmin;
            }
        }
    }
}

void Structure::add_residue(std::shared_ptr<io::IOResidue> residue) {
    if (residue) {
        residues_.push_back(residue);
    }
}

void Structure::add_atom(std::shared_ptr<io::PDBAtom> atom) {
    if (atom) {
        atoms_.push_back(atom);
    }
}

std::vector<std::array<double, 3>> Structure::get_coordinates() const {
    std::vector<std::array<double, 3>> coords;
    coords.reserve(atoms_.size());
    
    for (const auto& atom : atoms_) {
        if (atom) {
            coords.push_back({atom->x, atom->y, atom->z});
        }
    }
    
    return coords;
}

std::array<std::array<double, 3>, 3> Structure::get_box_vectors() const {
    std::array<std::array<double, 3>, 3> vectors{};
    
    if (box_) {
        const auto& [a, b, c, alpha, beta, gamma] = box_.value();
        
        // Convert angles from degrees to radians
        const double alpha_rad = alpha * M_PI / 180.0;
        const double beta_rad = beta * M_PI / 180.0;
        const double gamma_rad = gamma * M_PI / 180.0;
        
        // Calculate box vectors following the triclinic box convention
        vectors[0] = {a, 0.0, 0.0};
        vectors[1] = {b * cos(gamma_rad), b * sin(gamma_rad), 0.0};
        
        const double cx = c * cos(beta_rad);
        const double cy = c * (cos(alpha_rad) - cos(beta_rad) * cos(gamma_rad)) / sin(gamma_rad);
        const double cz = sqrt(c * c - cx * cx - cy * cy);
        vectors[2] = {cx, cy, cz};
    }
    
    return vectors;
}

std::unordered_map<std::string, double> Structure::get_energy_components() const {
    std::unordered_map<std::string, double> components;
    
    // Initialize all required components
    components["bond"] = 0.0;
    components["angle"] = 0.0;
    components["dihedral"] = 0.0;
    components["improper"] = 0.0;
    components["vdw"] = 0.0;
    components["electrostatic"] = 0.0;
    
    // Calculate energies if force field is available
    if (forcefield_) {
        // TODO: Implement actual energy calculations
        // For now, return placeholder values
        components["bond"] = 100.0;  // Example value
        components["angle"] = 200.0;  // Example value
        components["dihedral"] = 150.0;  // Example value
        components["improper"] = 50.0;  // Example value
        components["vdw"] = 300.0;  // Example value
        components["electrostatic"] = 400.0;  // Example value
    }
    
    // Calculate total energy
    double total = 0.0;
    for (const auto& [component, energy] : components) {
        if (component != "total") {
            total += energy;
        }
    }
    components["total"] = total;
    
    return components;
}

std::vector<std::tuple<size_t, double, double, double>> Structure::get_atom_energy_contributions() const {
    std::vector<std::tuple<size_t, double, double, double>> atom_energies;
    atom_energies.reserve(atoms_.size());
    
    if (!forcefield_) {
        return atom_energies;
    }
    
    // Calculate per-atom energy contributions
    for (size_t i = 0; i < atoms_.size(); ++i) {
        double vdw_energy = 0.0;
        double elec_energy = 0.0;
        
        const auto& atom1 = atoms_[i];
        if (!atom1) continue;
        
        for (size_t j = 0; j < atoms_.size(); ++j) {
            if (i == j) continue;
            
            const auto& atom2 = atoms_[j];
            if (!atom2) continue;
            
            // Calculate distance
            const double dx = atom2->x - atom1->x;
            const double dy = atom2->y - atom1->y;
            const double dz = atom2->z - atom1->z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r = sqrt(r2);
            
            // VDW energy contribution
            if (!std::isnan(atom1->forcefield_epsilon) && !std::isnan(atom2->forcefield_epsilon) &&
                !std::isnan(atom1->forcefield_rmin) && !std::isnan(atom2->forcefield_rmin)) {
                const double eps = sqrt(atom1->forcefield_epsilon * atom2->forcefield_epsilon);
                const double rmin = atom1->forcefield_rmin + atom2->forcefield_rmin;
                const double ratio = rmin / r;
                const double ratio6 = ratio * ratio * ratio * ratio * ratio * ratio;
                vdw_energy += 0.5 * eps * (ratio6 * ratio6 - 2.0 * ratio6);  // Half because of double counting
            }
            
            // Electrostatic energy contribution
            if (!std::isnan(atom1->topo_charge) && !std::isnan(atom2->topo_charge)) {
                const double k_coulomb = 332.0716;  // kcal/mol/Å/e^2
                elec_energy += 0.5 * k_coulomb * atom1->topo_charge * atom2->topo_charge / r;  // Half because of double counting
            }
        }
        
        atom_energies.emplace_back(i, vdw_energy, elec_energy, vdw_energy + elec_energy);
    }
    
    return atom_energies;
}

} // namespace core
} // namespace pygcmc 