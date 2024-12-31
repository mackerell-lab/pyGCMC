// modules/core/include/pygcmc/core/structure.hpp

#pragma once

#include "pygcmc/core/io/pdb_parser.hpp"
#include <memory>
#include <vector>

namespace pygcmc {
namespace core {

class ForceField;

class Structure {
public:
    Structure();
    ~Structure();

    // Apply force field parameters to atoms
    void apply_forcefield(std::shared_ptr<ForceField> ff);

    // Add residue and its atoms to the structure
    void add_residue(std::shared_ptr<io::IOResidue> residue);

    // Add a single atom to the structure
    void add_atom(std::shared_ptr<io::PDBAtom> atom);

    // Get number of atoms
    size_t get_num_atoms() const;

    // Access residues and atoms (const references to prevent modification)
    const std::vector<std::shared_ptr<io::IOResidue>>& residues() const;
    const std::vector<std::shared_ptr<io::PDBAtom>>& atoms() const;

private:
    std::vector<std::shared_ptr<io::IOResidue>> residues_;
    std::vector<std::shared_ptr<io::PDBAtom>> atoms_;
};

} // namespace core
} // namespace pygcmc 