#pragma once

#include <pygcmc/energy/energy.hpp>

namespace pygcmc {

class LennardJonesEnergy : public Energy {
public:
    LennardJonesEnergy() = default;
    ~LennardJonesEnergy() override = default;

    double calculate(const Project& project) override;
    double calculate_atom_energy(const Project& project, const ProjectAtom& atom) override;

private:
    // Helper function to calculate LJ potential between two atoms
    double calculate_pair_energy(const ProjectAtom& atom1, const ProjectAtom& atom2) const;
};

} // namespace pygcmc 