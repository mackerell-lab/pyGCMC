// modules/core/include/pygcmc/core/structure.hpp

#pragma once

#include <memory>
#include <vector>
#include <optional>
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
    void set_box(const std::optional<std::vector<double>>& box) { box_ = box; }
    std::optional<std::vector<double>> get_box() const { return box_; }

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

private:
    std::vector<std::shared_ptr<io::IOResidue>> residues_;
    std::vector<std::shared_ptr<io::PDBAtom>> atoms_;
    std::optional<std::vector<double>> box_;  // Box dimensions (a, b, c, alpha, beta, gamma)
};

} // namespace core
} // namespace pygcmc 