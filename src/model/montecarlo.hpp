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
 * - Atom:        Atomic coordinates (nm), type, charge (e) etc.
 * - Residue:     Group of atoms with start/count, for insertion/deletion/movement (swap-and-pop)
 * - GCMCSystem:  Contains pointers to (Atoms/Residues/ForceField/Info etc.)
 *                Provides interfaces for allocation, deallocation, download, upload
 *
 * Units used in this system follow GROMACS MD units, with conversions from PDB/CHARMM:
 * - Length: nanometers (nm) [converted from PDB/CHARMM Å]
 * - Energy: kilojoules per mole (kJ/mol) [converted from CHARMM kcal/mol]
 * - Charge: electron charge (e) [same as CHARMM partial charges]
 * - Time: picoseconds (ps)
 * - Temperature: Kelvin (K)
 * - Concentration: moles per liter (mol/L)
 * 
 * Unit Conversions:
 * - Distance: 1 Å = 0.1 nm
 * - Energy: 1 kcal/mol = 4.184 kJ/mol
 */

namespace pygcmc {
namespace model {

/**
 * @brief Type mapping system for atom and residue types
 * 
 * Provides bidirectional mapping between string type names and integer indices.
 * Used for both atom types and residue types in the system.
 */
struct TypeMaps {
    std::vector<std::string> atomTypes;  ///< Index -> Type string mapping
    std::unordered_map<std::string, int> atomTypeIndices;  ///< Type string -> Index mapping
    
    /**
     * @brief Get or add a type to the mapping
     * @param type Type string to look up or add
     * @return Index of the type in the mapping
     */
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
    
    /**
     * @brief Get type name from index
     * @param index Index to look up
     * @return Type string, empty if index invalid
     */
    std::string getTypeName(int index) const {
        if (index >= 0 && static_cast<size_t>(index) < atomTypes.size()) {
            return atomTypes[index];
        }
        return "";
    }
};

/**
 * @brief Global Monte Carlo simulation parameters
 * 
 * Contains all global parameters needed for the MC simulation,
 * including system dimensions, thermodynamic conditions, and statistics.
 */
struct MCInfo {
    int    mcSteps{0};       ///< Monte Carlo steps (dimensionless)
    float  box[3]{-1.0f};    ///< Box dimensions (nm, converted from PDB CRYST1 record in Å)
    float  cutoff{1.5f};     ///< Cutoff distance for non-bonded interactions (nm)
    
    /// @brief Inverse temperature beta = 1/(kB*T) [mol/kJ]
    /// @note Used for Metropolis criterion with energies in kJ/mol
    float  beta{0.0f};       

    // Reserved max capacity
    int    maxResidues;      ///< Maximum number of residues
    int    maxAtoms;         ///< Maximum number of atoms
    int    maxTypes;         ///< Maximum number of atom types

    // Global parameters
    float  volume;           ///< System volume (nm³)
    uint64_t seed;          ///< Random seed

    /**
     * @brief Statistics for Monte Carlo moves
     * 
     * Tracks acceptance rates for different types of MC moves
     */
    struct Statistics {
        int totalMoves{0};           ///< Total number of moves attempted
        int acceptedMoves{0};        ///< Number of accepted moves
        int insertionAttempts{0};    ///< Number of insertion attempts
        int acceptedInsertions{0};   ///< Number of accepted insertions
        int deletionAttempts{0};     ///< Number of deletion attempts
        int acceptedDeletions{0};    ///< Number of accepted deletions
    } stats;

    // Constants (GROMACS MD units)
    /// @brief Boltzmann constant [kJ/(mol·K)]
    /// @note Consistent with energy unit conversion from CHARMM (kcal/mol -> kJ/mol)
    static constexpr float BOLTZMANN = 0.00831446f;  
    
    /// @brief Energy conversion factor from CHARMM to internal units
    static constexpr float KCAL_TO_KJ = 4.184f;      
    static constexpr float KJ_TO_KCAL = 0.239f;      
    static constexpr float MOLES_TO_MOLECULES = 0.0006023f;  
    static constexpr float MOLECULES_TO_MOLES = 1660.539f;   

    /**
     * @brief Set temperature and calculate beta
     * @param temperature Temperature in Kelvin
     */
    void setTemperature(float temperature) {
        beta = 1.0f / (BOLTZMANN * temperature);
    }
};

/**
 * @brief Force field parameters for Monte Carlo simulation
 * 
 * CHARMM force field conversion:
 * - Original CHARMM format uses Rmin/2 (Å) and epsilon (kcal/mol)
 * - Internal format uses sigma (nm) and epsilon (kJ/mol)
 * 
 * Conversions:
 * 1. Distance: 
 *    - CHARMM provides Rmin/2 in Å
 *    - Convert to nm: multiply by 0.1
 * 
 * 2. Energy:
 *    - CHARMM provides epsilon in kcal/mol
 *    - Convert to kJ/mol: multiply by 4.184
 * 
 * Lennard-Jones potential form:
 * - E = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
 * 
 * Storage layout:
 * - Arrays are 1D but represent 2D interaction matrices
 * - Size is (numTotalTypes * numTotalTypes)
 * - For type i and j: index = i * numTotalTypes + j
 */
struct MCForceField {
    /// @brief Total number of atom types in system
    int numTotalTypes;    

    /// @brief Number of atom types in movement molecules
    /// @note A single movement molecule may contain multiple atom types
    /// @deprecated This field will be used in future optimization
    int numMovementTypes;   

    /// @brief LJ sigma parameters [nm] converted from CHARMM Rmin/2 [Å]
    /// @note sigma = (Rmin/2) * 2 * 2^(-1/6) * 0.1
    /// @note Stores all type pairs as (total_type, total_type) matrix
    /// @note Size is numTotalTypes * numTotalTypes
    std::vector<float> ljSigma;   

    /// @brief LJ epsilon parameters [kJ/mol] converted from CHARMM [kcal/mol]
    /// @note epsilon_kj = epsilon_kcal * 4.184
    /// @note Stores all type pairs as (total_type, total_type) matrix
    /// @note Size is numTotalTypes * numTotalTypes
    std::vector<float> ljEps;     
};

/**
 * @brief Basic atomic properties for Monte Carlo simulation
 * 
 * Represents a single atom in the system with its position,
 * charge, and type information.
 */
struct MCAtom {
    float x, y, z;      ///< Position (nm, converted from Å in PDB)
    float charge;       ///< Charge (e, CHARMM partial charges)
    int   type;        ///< Type index mapping to CHARMM atom types
};

/**
 * @brief Molecular unit for GCMC simulation
 * 
 * Represents a group of atoms that move together in the simulation.
 * Used for insertion, deletion, and movement operations.
 */
struct MCResidue {
    // Basic properties
    int   atomStart;    ///< Starting index in global atom array
    int   atomCount;    ///< Number of atoms in this residue
    bool  active;       ///< Whether this residue is currently in use
    bool  fixed;        ///< Whether this residue can be moved
    float center[3];    ///< Geometric center (nm, derived from atom positions)

    // Energy components (converted from CHARMM kcal/mol to kJ/mol)
    float energy_vdw;   ///< Lennard-Jones energy (kJ/mol)
    float energy_elec;  ///< Coulomb energy (kJ/mol)
        
    // GCMC parameters
    /// @brief Target concentration for GCMC insertion/deletion (mol/L)
    /// @note Used to determine acceptance probability
    float concentration; 

    /// @brief Excess chemical potential for GCMC moves (kJ/mol)
    /// @note Should be consistent with CHARMM force field energies after conversion
    float chemPot;      

    int   type;         ///< Residue type index mapping to CHARMM residue names
    float radius;       ///< Approximate radius (nm)
};

/**
 * @brief Information about movement residues in the system
 * 
 * Tracks the location and count of residues that can be moved,
 * inserted, or deleted during the simulation.
 */
struct MCMovementResidueInfo {
    int startIndex;         ///< Starting index in global residue array
    int activeCount;        ///< Number of active movement residues
    int totalCount;         ///< Total number of movement residues (active + inactive)
    std::string resName;    ///< Residue name for this movement group
};

/**
 * @brief Current state of the Monte Carlo system
 * 
 * Contains all information about the current state of the system,
 * including atoms, residues, type mappings, and force field parameters.
 */
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
