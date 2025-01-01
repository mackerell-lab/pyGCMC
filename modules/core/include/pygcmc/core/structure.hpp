// modules/core/include/pygcmc/core/structure.hpp

#pragma once

#include <memory>
#include <vector>
#include <optional>
#include <array>
#include <unordered_map>
#include <tuple>
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/forcefield.hpp"

namespace pygcmc {
namespace core {

class ForceField;

class Structure {
public:
    Structure() = default;
    ~Structure() = default;

    // Box dimensions methods
    void set_box(const std::optional<std::array<double, 6>>& box) { box_ = box; }
    std::optional<std::array<double, 6>> get_box() const { return box_; }

    // Apply force field parameters to all atoms
    void apply_forcefield(const std::shared_ptr<ForceField>& forcefield);

    // Add a residue to the structure
    void add_residue(std::shared_ptr<io::IOResidue> residue);

    // Add an atom to the structure
    void add_atom(std::shared_ptr<io::PDBAtom> atom);

    // Get number of atoms
    size_t get_num_atoms() const { return atoms_.size(); }

    // Get residues
    const std::vector<std::shared_ptr<io::IOResidue>>& residues() const { return residues_; }

    // Get atoms
    const std::vector<std::shared_ptr<io::PDBAtom>>& atoms() const { return atoms_; }

    // New data access methods
    std::vector<std::array<double, 3>> get_coordinates() const;
    std::array<std::array<double, 3>, 3> get_box_vectors() const;
    
    // Energy calculation methods
    std::unordered_map<std::string, double> get_energy_components() const;
    std::vector<std::tuple<size_t, double, double, double>> get_atom_energy_contributions() const;

private:
    std::vector<std::shared_ptr<io::IOResidue>> residues_;
    std::vector<std::shared_ptr<io::PDBAtom>> atoms_;
    std::optional<std::array<double, 6>> box_;  // Box dimensions (a, b, c, alpha, beta, gamma)
    std::shared_ptr<ForceField> forcefield_;
};

} // namespace core
} // namespace pygcmc 