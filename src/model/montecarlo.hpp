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
struct MCInfo {
    int    mcSteps{0};       ///< Monte Carlo steps
    float  box[3]{-1.0f};        ///< Box dimensions (Å)
    float  cutoff{15.0f};       ///< Cutoff distance for non-bonded interactions (Å)
    float  beta{0.0f};          ///< 1/(kB*T) (mol/kcal) - will be set based on temperature
    
    // Reserved max capacity
    int    maxResidues;   
    int    maxAtoms;
    int    maxTypes;      ///< Maximum number of atom types

    // Global parameters
    float  volume;        ///< System volume
    uint64_t seed;       ///< Random seed

    // Statistics
    struct Statistics {
        int totalMoves{0};           ///< Total number of moves attempted
        int acceptedMoves{0};        ///< Number of accepted moves
        int insertionAttempts{0};    ///< Number of insertion attempts
        int acceptedInsertions{0};   ///< Number of accepted insertions
        int deletionAttempts{0};     ///< Number of deletion attempts
        int acceptedDeletions{0};    ///< Number of accepted deletions
    } stats;

    // Constants (CHARMM units)
    static constexpr float BOLTZMANN = 0.0019881f;  ///< Boltzmann constant (kcal/mol/K)
    static constexpr float KCAL_TO_KJ = 4.184f;     ///< Convert kcal/mol to kJ/mol
    static constexpr float KJ_TO_KCAL = 0.239f;     ///< Convert kJ/mol to kcal/mol
    static constexpr float MOLES_TO_MOLECULES = 0.0006023f;  ///< Convert mol/L to molecules/Å³
    static constexpr float MOLECULES_TO_MOLES = 1660.539f;   ///< Convert molecules/Å³ to mol/L

    // Set temperature and calculate beta
    void setTemperature(float temperature) {  // temperature in Kelvin
        beta = 1.0f / (BOLTZMANN * temperature);
    }
};

// ------------------------------------------------------------
// 2) ForceField: Force field parameters
// ------------------------------------------------------------
struct MCForceField {
    int maxTypes;  ///< Actual number of atom types in use

    // Lennard-Jones parameters
    std::vector<float> ljSigma;   ///< [maxTypes] sigma parameters
    std::vector<float> ljEps;     ///< [maxTypes] epsilon parameters
};

// ------------------------------------------------------------
// 3) Atom: Basic atomic properties
// ------------------------------------------------------------
struct MCAtom {
    float x, y, z;      ///< Position
    float charge;       ///< Charge
    int   type;        ///< Type index
};

// ------------------------------------------------------------
// 4) Residue: Molecular unit for GCMC
// ------------------------------------------------------------
struct MCResidue {
    // Basic properties
    int   atomStart;    ///< Starting index in global atom array
    int   atomCount;    ///< Number of atoms
    bool  active;       ///< Whether in use
    bool  fixed;        ///< Whether fixed  
    float center[3];    ///< Geometric center
    
    // GCMC parameters
    float concentration; ///< Target concentration
    float chemPot;      ///< Chemical potential
    int   type;         ///< Residue type
    float radius;       ///< Approximate radius
};

// ------------------------------------------------------------
// 4.5) Movement Residue Info: Track GCMC movement residues
// ------------------------------------------------------------
struct MCMovementResidueInfo {
    int startIndex;         ///< Starting index of movement residues in global residue array
    int activeCount;        ///< Number of active movement residues
    int totalCount;         ///< Total number of movement residues (active + inactive)
    std::string resName;    ///< Residue name for this movement group
};

// ------------------------------------------------------------
// 5) System State: Current state of the MC system
// ------------------------------------------------------------
struct MCState {
    // Arrays
    std::vector<MCAtom>    atoms;      ///< Global atom array
    std::vector<MCResidue> residues;   ///< Global residue array

    // Type mappings
    TypeMaps residueTypes;  ///< Residue type mappings
    TypeMaps atomTypes;     ///< Atom type mappings

    // Active counts for swap-and-pop management
    int activeAtomCount;     ///< Current active atom count
    int activeResidueCount;  ///< Current active residue count

    // Movement residue tracking
    std::vector<MCMovementResidueInfo> movementResidues;  ///< Info for each type of movement residue

    // Parameters
    MCInfo   info;        ///< Global MC parameters
    MCForceField forcefield;  ///< Force field parameters

    // Constructor to initialize counts
    MCState() : activeAtomCount(0), activeResidueCount(0) {}
};

} // namespace gcmc
