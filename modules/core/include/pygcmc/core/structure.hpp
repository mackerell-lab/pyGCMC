// modules/core/include/pygcmc/core/structure.hpp

#pragma once

#include <memory>
#include <vector>
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/forcefield.hpp"

namespace pygcmc {
namespace core {

class ForceField;

class Structure {
public:
    Structure();
    ~Structure();

    // Apply force field parameters to all atoms
    void apply_forcefield(const std::shared_ptr<ForceField>& forcefield);

    // Add a residue to the structure
    void add_residue(std::shared_ptr<io::IOResidue> residue);

    // Add an atom to the structure
    void add_atom(std::shared_ptr<io::PDBAtom> atom);

    // Get number of atoms
    size_t get_num_atoms() const;

    // Get residues
    const std::vector<std::shared_ptr<io::IOResidue>>& residues() const;

    // Get atoms
    const std::vector<std::shared_ptr<io::PDBAtom>>& atoms() const;

private:
    std::vector<std::shared_ptr<io::IOResidue>> residues_;
    std::vector<std::shared_ptr<io::PDBAtom>> atoms_;
};

} // namespace core
} // namespace pygcmc 