#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

/**
 * @file   montecarlo.hpp
 * @brief  Core data structures for GCMC simulation: GCMCInfo, ForceField, Atom, Residue + GCMCSystem manager
 *
 * - GCMCInfo:    Global MC parameters (temperature, box, max molecules/atoms etc.)
 * - ForceField:  Force field parameters (e.g. Lennard-Jones), extensible for bond/angle/dihedral
 * - Atom:        Atomic coordinates, type, charge etc.
 * - Residue:     Group of atoms with start/count, for insertion/deletion/movement (swap-and-pop)
 * - GCMCSystem:  Contains pointers to (Atoms/Residues/ForceField/Info etc.)
 *                Provides interfaces for allocation, deallocation, download, upload
 */

namespace gcmc {

// ------------------------------------------------------------
// 0) TypeMaps: Atom type mapping system
// ------------------------------------------------------------
struct TypeMaps {
    std::vector<std::string> atomTypes;  // Index -> Type string mapping
    std::unordered_map<std::string, int> atomTypeIndices;  // Type string -> Index mapping
    
    int getOrAddType(const std::string& type) {
        auto it = atomTypeIndices.find(type);
        if (it != atomTypeIndices.end()) {
            return it->second;
        }
        int newIndex = atomTypes.size();
        atomTypes.push_back(type);
        atomTypeIndices[type] = newIndex;
        return newIndex;
    }
    
    std::string getTypeName(int index) const {
        if (index >= 0 && static_cast<size_t>(index) < atomTypes.size()) {
            return atomTypes[index];
        }
        return "";
    }
};

// ------------------------------------------------------------
// 1) GCMCInfo: Global MC parameters
// ------------------------------------------------------------
struct GCMCInfo {
    int    mcSteps;       ///< Monte Carlo steps
    float  box[3];        ///< Box dimensions (x, y, z)
    float  cutoff;        ///< Cutoff distance (Å)
    float  beta;          ///< 1/(kB*T)
    
    // Reserved max capacity
    int    maxResidues;   
    int    maxAtoms;
    int    maxTypes;      ///< Maximum number of atom types

    // Global parameters
    float  volume;        ///< System volume
    uint64_t seed;       ///< Random seed

    // Statistics
    struct Statistics {
        int totalMoves;           ///< Total number of moves attempted
        int acceptedMoves;        ///< Number of accepted moves
        int insertionAttempts;    ///< Number of insertion attempts
        int acceptedInsertions;   ///< Number of accepted insertions
        int deletionAttempts;     ///< Number of deletion attempts
        int acceptedDeletions;    ///< Number of accepted deletions
    } stats;
};

// ------------------------------------------------------------
// 2) ForceField: Force field parameters
// ------------------------------------------------------------
struct ForceField {
    int maxTypes;  ///< Actual number of atom types in use

    // Lennard-Jones parameters
    std::vector<float> ljSigma;   ///< [maxTypes] sigma parameters
    std::vector<float> ljEps;     ///< [maxTypes] epsilon parameters
};

// ------------------------------------------------------------
// 3) Atom: Basic atomic properties
// ------------------------------------------------------------
struct Atom {
    float x, y, z;      ///< Position
    float charge;       ///< Charge
    int   type;        ///< Type index
};

// ------------------------------------------------------------
// 4) Residue: Molecular unit for GCMC
// ------------------------------------------------------------
struct Residue {
    // Basic properties
    int   atomStart;    ///< Starting index in global atom array
    int   atomCount;    ///< Number of atoms
    bool  active;       ///< Whether in use
    float com[3];       ///< Center of mass
    
    // GCMC parameters
    float concentration; ///< Target concentration
    float chemPot;      ///< Chemical potential
    int   type;         ///< Residue type
    float radius;       ///< Approximate radius
};

// ------------------------------------------------------------
// 5) System State: Current state of the MC system
// ------------------------------------------------------------
struct SystemState {
    // Arrays
    std::vector<Atom>    atoms;      ///< Global atom array
    std::vector<Residue> residues;   ///< Global residue array

    // Active counts for swap-and-pop management
    int activeAtomCount;     ///< Current active atom count
    int activeResidueCount;  ///< Current active residue count

    // Parameters
    GCMCInfo   info;        ///< Global MC parameters
    ForceField forcefield;  ///< Force field parameters

    // Constructor to initialize counts
    SystemState() : activeAtomCount(0), activeResidueCount(0) {}
};

} // namespace gcmc
