// modules/core/src/structure.cpp

#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace pygcmc {
namespace core {

Structure::Structure() {}

Structure::~Structure() {
    // Clear vectors in reverse order of dependency
    atoms_.clear();
    residues_.clear();
}

void Structure::apply_forcefield(const std::shared_ptr<ForceField>& forcefield) {
    if (!forcefield) {
        std::cerr << "Error: Null force field pointer provided to apply_forcefield" << std::endl;
        return;
    }

    std::cout << "\nApplying force field parameters to atoms:" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::left 
              << std::setw(10) << "Residue"
              << std::setw(6) << "Seq"
              << std::setw(8) << "Name"
              << std::setw(10) << "TopoType"
              << std::setw(15) << "Old Epsilon"
              << std::setw(15) << "New Epsilon"
              << std::setw(15) << "New Rmin" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    for (auto& atom_ptr : atoms_) {
        if (!atom_ptr) continue;
        
        const auto& topo_type = atom_ptr->topo_type;
        const auto old_epsilon = atom_ptr->forcefield_epsilon;
        
        auto it = forcefield->nonbonded_params().find(topo_type);
        if (it != forcefield->nonbonded_params().end()) {
            atom_ptr->forcefield_epsilon = it->second.epsilon;
            atom_ptr->forcefield_rmin = it->second.rmin;
            
            std::cout << std::left 
                      << std::setw(10) << atom_ptr->residue
                      << std::setw(6) << atom_ptr->sequence
                      << std::setw(8) << atom_ptr->name
                      << std::setw(10) << topo_type
                      << std::setw(15) << old_epsilon
                      << std::setw(15) << atom_ptr->forcefield_epsilon
                      << std::setw(15) << atom_ptr->forcefield_rmin << std::endl;
        } else {
            std::cerr << "Warning: No force field parameters found for topo_type '" 
                      << topo_type << "' in atom " << atom_ptr->name 
                      << " of residue " << atom_ptr->residue 
                      << " " << atom_ptr->sequence << std::endl;
        }
    }
    std::cout << std::endl;
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
        
        // Add residue to residues list first
        residues_.push_back(residue);
        
        // Add all atoms from this residue's atom_ptrs
        for (const auto& atom_ptr : residue->atom_ptrs) {
            if (!atom_ptr || !atom_ptr->is_valid()) {
                // Rollback on error
                residues_.pop_back();
                throw std::invalid_argument("Invalid atom in residue");
            }
            atoms_.push_back(atom_ptr);  // Share ownership of the atom
        }
    } catch (...) {
        throw;  // Re-throw after cleanup
    }
}

void Structure::add_atom(std::shared_ptr<io::PDBAtom> atom) {
    if (!atom || !atom->is_valid()) {
        throw std::invalid_argument("Invalid atom");
    }
    atoms_.push_back(atom);  // Share ownership of the atom
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