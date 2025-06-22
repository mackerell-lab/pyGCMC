#pragma once

#ifndef PYGCMC_MODEL_UTILS_FACTORY_HPP
#define PYGCMC_MODEL_UTILS_FACTORY_HPP

#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../molecule/MolecularMain.hpp"
#include <memory>
#include <string>

namespace pygcmc {
namespace model {
namespace factory {

// Using declarations for cleaner code
using AtomPtr = std::shared_ptr<atom::Atom>;
using ResiduePtr = std::shared_ptr<residue::Residue>;
using MolecularPtr = std::shared_ptr<molecule::Molecular>;

/**
 * @brief Create a water molecule
 * @param residue_id The residue ID for the water molecule
 * @param segment The segment name (default: "SOLV")
 * @return Shared pointer to the created water residue
 */
inline ResiduePtr create_water_molecule(int residue_id, const std::string& segment = "SOLV") {
    auto water = std::make_shared<residue::Residue>("SOL", residue_id, segment);
    
    // Add atoms: O, H1, H2 with correct constructor parameters
    // AtomCore(int bynu, type, resname, ires, segid, iseg, x, y, z, wmain, mass, charge, chem, hetatm)
    auto oxygen = std::make_shared<atom::Atom>(1, "OT", "SOL", residue_id, segment, 0,
                                               0.0, 0.0, 0.0, 1.0, 15.999, -0.834, "O", false);
    auto h1 = std::make_shared<atom::Atom>(2, "HT", "SOL", residue_id, segment, 0,
                                           0.757, 0.586, 0.0, 1.0, 1.008, 0.417, "H", false);
    auto h2 = std::make_shared<atom::Atom>(3, "HT", "SOL", residue_id, segment, 0,
                                           -0.757, 0.586, 0.0, 1.0, 1.008, 0.417, "H", false);
    
    // Add atoms to residue
    water->add_atom(oxygen);
    water->add_atom(h1);
    water->add_atom(h2);
    
    return water;
}

/**
 * @brief Create a simple protein residue (Alanine)
 * @param residue_id The residue ID for the alanine residue
 * @param segment The segment name (default: "PROT")
 * @return Shared pointer to the created alanine residue
 */
inline ResiduePtr create_alanine_residue(int residue_id, const std::string& segment = "PROT") {
    auto ala = std::make_shared<residue::Residue>("ALA", residue_id, segment);
    
    // Add backbone atoms with correct constructor parameters
    auto n = std::make_shared<atom::Atom>(1, "NH1", "ALA", residue_id, segment, 0,
                                          0.0, 0.0, 0.0, 1.0, 14.007, -0.47, "N", false);
    auto ca = std::make_shared<atom::Atom>(2, "CT1", "ALA", residue_id, segment, 0,
                                           1.458, 0.0, 0.0, 1.0, 12.01, 0.07, "C", false);
    auto c = std::make_shared<atom::Atom>(3, "C", "ALA", residue_id, segment, 0,
                                          2.009, 1.421, 0.0, 1.0, 12.01, 0.51, "C", false);
    auto o = std::make_shared<atom::Atom>(4, "O", "ALA", residue_id, segment, 0,
                                          1.239, 2.364, 0.0, 1.0, 15.999, -0.51, "O", false);
    auto cb = std::make_shared<atom::Atom>(5, "CT3", "ALA", residue_id, segment, 0,
                                           2.196, -0.889, -1.07, 1.0, 12.01, -0.27, "C", false);
    
    // Add atoms to residue
    ala->add_atom(n);
    ala->add_atom(ca);
    ala->add_atom(c);
    ala->add_atom(o);
    ala->add_atom(cb);
    
    return ala;
}

/**
 * @brief Create a molecular system with given residues
 * @param residues Vector of residue pointers to add to the system
 * @return Shared pointer to the created molecular system
 */
inline MolecularPtr create_molecular_system(const std::vector<ResiduePtr>& residues) {
    auto molecular = std::make_shared<molecule::Molecular>();
    
    for (const auto& residue : residues) {
        if (residue) {
            molecular->add_residue(residue);
        }
    }
    
    return molecular;
}

/**
 * @brief Create a simple molecular system with water and protein residues
 * @param num_waters Number of water molecules to create
 * @param num_alanines Number of alanine residues to create
 * @return Shared pointer to the created molecular system
 */
inline MolecularPtr create_simple_system(int num_waters = 10, int num_alanines = 1) {
    std::vector<ResiduePtr> residues;
    
    // Add water molecules
    for (int i = 0; i < num_waters; ++i) {
        residues.push_back(create_water_molecule(i + 1, "SOLV"));
    }
    
    // Add alanine residues
    for (int i = 0; i < num_alanines; ++i) {
        residues.push_back(create_alanine_residue(num_waters + i + 1, "PROT"));
    }
    
    return create_molecular_system(residues);
}

} // namespace factory
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_UTILS_FACTORY_HPP 