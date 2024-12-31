// modules/core/src/structure.cpp

#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"
#include <memory>
#include <algorithm>

namespace pygcmc {
namespace core {

Structure::Structure() {}

Structure::~Structure() {
    // Clear vectors in reverse order of dependency
    atoms_.clear();
    residues_.clear();
}

void Structure::apply_forcefield(std::shared_ptr<ForceField> ff) {
    if (!ff) return;  // Skip if force field is null
    
    // Apply force field parameters to all atoms
    for (auto& atom : atoms_) {
        if (!atom) continue;  // Skip null pointers
        
        try {
            // Get nonbonded parameters for this atom type
            auto it = ff->get_nonbonded_params().find(atom->type);
            if (it != ff->get_nonbonded_params().end()) {
                atom->forcefield_epsilon = it->second.epsilon;
                atom->forcefield_rmin = it->second.rmin;
            }
        } catch (const std::exception& e) {
            // Log error but continue with other atoms
            continue;
        }
    }
}

void Structure::add_residue(std::shared_ptr<io::IOResidue> residue) {
    if (!residue) return;  // Skip null pointers
    
    try {
        // Validate residue before adding
        if (!residue->is_valid()) {
            throw std::invalid_argument("Invalid residue");
        }
        
        // Pre-calculate new size to check for overflow
        const size_t current_atoms = atoms_.size();
        const size_t new_atoms = residue->atom_ptrs.size();
        if (new_atoms > std::numeric_limits<size_t>::max() - current_atoms) {
            throw std::overflow_error("Too many atoms");
        }
        
        // Reserve space for new atoms
        atoms_.reserve(current_atoms + new_atoms);
        
        // Keep track of successfully added atoms for rollback
        const size_t original_size = atoms_.size();
        bool added_residue = false;
        
        try {
            // Add residue to residues list first
            residues_.push_back(residue);
            added_residue = true;
            
            // Add all atoms from this residue's atom_ptrs
            for (const auto& atom_ptr : residue->atom_ptrs) {
                if (!atom_ptr || !atom_ptr->is_valid()) {
                    throw std::invalid_argument("Invalid atom in residue");
                }
                // Create a new shared_ptr that shares ownership
                atoms_.push_back(atom_ptr);
            }
        } catch (...) {
            // Rollback on any error
            if (added_residue) {
                residues_.pop_back();
            }
            // Restore atoms_ to original state
            atoms_.resize(original_size);
            throw;  // Re-throw the exception
        }
    } catch (const std::exception& e) {
        throw;  // Re-throw after cleanup
    }
}

void Structure::add_atom(std::shared_ptr<io::PDBAtom> atom) {
    if (!atom) return;  // Skip null pointers
    if (!atom->is_valid()) {
        throw std::invalid_argument("Invalid atom");
    }
    atoms_.push_back(atom);
}

size_t Structure::get_num_atoms() const {
    return atoms_.size();
}

const std::vector<std::shared_ptr<io::IOResidue>>& Structure::residues() const {
    return residues_;
}

const std::vector<std::shared_ptr<io::PDBAtom>>& Structure::atoms() const {
    return atoms_;
}

} // namespace core
} // namespace pygcmc 