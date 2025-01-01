// modules/core/src/structure.cpp

#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"
#include <memory>
#include <vector>
#include <stdexcept>

namespace pygcmc {
namespace core {

void Structure::apply_forcefield(const std::shared_ptr<ForceField>& forcefield) {
    if (!forcefield) {
        throw std::runtime_error("Cannot apply null forcefield");
    }
    
    // Apply forcefield parameters to each atom
    for (auto& atom : atoms_) {
        if (atom) {
            // Apply parameters based on atom type
            // ... implementation ...
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

} // namespace core
} // namespace pygcmc 