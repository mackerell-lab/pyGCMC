// src/model/montecarlo.hpp

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

namespace pygcmc {
namespace model {

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
    int    mcSteps{0};       ///< Monte Carlo steps (dimensionless)
    float  box[3]{-1.0f};    ///< Box dimensions (Å)
    float  cutoff{15.0f};    ///< Cutoff distance for non-bonded interactions (Å)
    float  beta{0.0f};       ///< 1/(kB*T) (mol/kcal) - inverse temperature
    
    // Reserved max capacity
    int    maxResidues;      ///< Maximum number of residues (dimensionless)
    int    maxAtoms;         ///< Maximum number of atoms (dimensionless)
    int    maxTypes;         ///< Maximum number of atom types (dimensionless)

    // Global parameters
    float  volume;           ///< System volume (Å³)
    uint64_t seed;          ///< Random seed (dimensionless)

    // Statistics
    struct Statistics {
        int totalMoves{0};           ///< Total number of moves attempted (dimensionless)
        int acceptedMoves{0};        ///< Number of accepted moves (dimensionless)
        int insertionAttempts{0};    ///< Number of insertion attempts (dimensionless)
        int acceptedInsertions{0};   ///< Number of accepted insertions (dimensionless)
        int deletionAttempts{0};     ///< Number of deletion attempts (dimensionless)
        int acceptedDeletions{0};    ///< Number of accepted deletions (dimensionless)
    } stats;

    // Constants (CHARMM units)
    static constexpr float BOLTZMANN = 0.0019881f;  ///< Boltzmann constant (kcal/mol/K)
    static constexpr float KCAL_TO_KJ = 4.184f;     ///< Convert kcal/mol to kJ/mol
    static constexpr float KJ_TO_KCAL = 0.239f;     ///< Convert kJ/mol to kcal/mol
    static constexpr float MOLES_TO_MOLECULES = 0.0006023f;  ///< Convert mol/L to molecules/Å³
    static constexpr float MOLECULES_TO_MOLES = 1660.539f;   ///< Convert molecules/Å³ to mol/L

    // Set temperature and calculate beta
    void setTemperature(float temperature) {  // temperature in Kelvin
        // beta = 1/(kB*T) where kB is BOLTZMANN in kcal/mol/K and T is in K
        // This gives beta in mol/kcal units
        beta = 1.0f / (BOLTZMANN * temperature);
    }
};

// ------------------------------------------------------------
// 2) ForceField: Force field parameters
// ------------------------------------------------------------
/**
 * @brief Force field parameters for Monte Carlo simulation
 * 
 * The force field parameters are organized to optimize the energy calculation between
 * movement molecules and all other molecules in the system. Parameters are stored in
 * 1D arrays but represent a 2D matrix of interactions:
 *   - Rows: movement atom types (numMovementTypes)
 *   - Columns: all possible atom types in the system (maxTypes)
 * 
 * For a movement atom type i and any atom type j, parameters are accessed using:
 * index = i * maxTypes + j
 * 
 * The Lennard-Jones potential is defined as:
 *   V(r) = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
 * where:
 *   - eps: well depth (in kcal/mole), stored in ljEps
 *   - sigma: distance at which potential is zero (in Angstroms), stored in ljSigma
 *   - r: distance between atoms
 * 
 * Parameters are combined using Lorentz-Berthelot rules:
 *   - sigma: arithmetic mean (Lorentz)
 *      sigma_ij = (sigma_i + sigma_j) / 2
 *   - epsilon: geometric mean (Berthelot)
 *      eps_ij = sqrt(eps_i * eps_j)
 * 
 * For NBFIX pairs, specific eps values are used directly, but sigma still uses
 * arithmetic mean.
 */
struct MCForceField {
    int maxTypes;           ///< Total number of atom types in the system
    int numMovementTypes;   ///< Number of atom types that belong to movement molecules

    // Arrays store parameters for movement types interacting with all types
    // Size: numMovementTypes * maxTypes
    std::vector<float> ljSigma;   ///< Combined sigma values [Å] for each type pair
    std::vector<float> ljEps;     ///< Combined epsilon values [kcal/mole] for each type pair
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
    float center[3];    ///< Geometric center (Å): arithmetic mean of all atom coordinates in this residue

    // Energy components
    float energy_vdw;   ///< Lennard-Jones energy (kcal/mole)
    float energy_elec;  ///< Coulomb energy (kcal/mole)
        
    // GCMC parameters
    float concentration; ///< Target concentration (mol/L)
    float chemPot;      ///< Chemical potential (kcal/mole)
    int   type;         ///< Residue type index in residueTypes map
    float radius;       ///< Approximate radius (Å)
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

    // Movement molecule info
    std::vector<MCMovementResidueInfo> movementResidues;  ///< Info for movement residues
    std::vector<int> movementAtomTypes;  ///< Atom types belonging to movement molecules
    int numMovementAtomTypes;  ///< Number of atom types from movement molecules

    // Active counts for swap-and-pop management
    int activeAtomCount;     ///< Current active atom count
    int activeResidueCount;  ///< Current active residue count

    // System info and force field
    MCInfo info;            ///< System information
    MCForceField forcefield;  ///< Force field parameters

    // Constructor to initialize counts
    MCState() : activeAtomCount(0), activeResidueCount(0) {}
};

} // namespace model
} // namespace pygcmc
