#pragma once

#include <pygcmc/energy/energy.hpp>

namespace pygcmc {

class ElectrostaticEnergy : public Energy {
public:
    ElectrostaticEnergy() = default;
    ~ElectrostaticEnergy() override = default;

    double calculate(const Project& project) override;
    double calculate_atom_energy(const Project& project, const ProjectAtom& atom) override;

private:
    // Helper function to calculate Coulomb potential between two atoms
    double calculate_pair_energy(const ProjectAtom& atom1, const ProjectAtom& atom2) const;
    
    // Coulomb constant in kcal⋅mol⁻¹⋅Å⋅e⁻²
    static constexpr double k_coulomb = 332.0716;
};

} // namespace pygcmc 