#pragma once

#ifndef PYGCMC_MODEL_MONTECARLO_STATE_CORE_HPP
#define PYGCMC_MODEL_MONTECARLO_STATE_CORE_HPP

#include "../common/ModelInterface.hpp"
#include "../common/ModelConstants.hpp"
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <array>

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Type mapping system for atom and residue types
 * Provides bidirectional mapping between string type names and integer indices
 */
struct TypeMaps {
    std::vector<std::string> atomTypes;  ///< Index -> Type string mapping
    std::unordered_map<std::string, int> atomTypeIndices;  ///< Type string -> Index mapping
    
    /**
     * @brief Get or add a type to the mapping
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
     */
    std::string getTypeName(int index) const {
        if (index >= 0 && static_cast<size_t>(index) < atomTypes.size()) {
            return atomTypes[index];
        }
        return "";
    }

    /**
     * @brief Get number of types
     */
    size_t size() const { return atomTypes.size(); }

    /**
     * @brief Clear all mappings
     */
    void clear() {
        atomTypes.clear();
        atomTypeIndices.clear();
    }
};

/**
 * @brief Global Monte Carlo simulation parameters
 * Contains all global parameters needed for the MC simulation
 */
struct MCInfo {
    int    mcSteps{0};       ///< Monte Carlo steps (dimensionless)
    float  box[3]{-1.0f};    ///< Box dimensions (nm)
    float  cutoff{1.5f};     ///< Cutoff distance for non-bonded interactions (nm)
    
    /// @brief Inverse temperature beta = 1/(kB*T) [mol/kJ]
    float  beta{0.0f};       

    // Reserved max capacity
    int    maxResidues;      ///< Maximum number of residues
    int    maxAtoms;         ///< Maximum number of atoms
    int    maxTypes;         ///< Maximum number of atom types

    // Global parameters
    float  volume;           ///< System volume (nm³)
    uint64_t seed;          ///< Random seed

    // CHARMM-style switching function parameters
    bool use_switching{false};  ///< Whether to use switching function
    float r_on{1.0f};           ///< Inner cutoff radius (nm)
    float r_off{1.2f};          ///< Outer cutoff radius (nm)

    /**
     * @brief Statistics for Monte Carlo moves
     */
    struct Statistics {
        int totalMoves{0};           ///< Total number of moves attempted
        int acceptedMoves{0};        ///< Number of accepted moves
        int insertionAttempts{0};    ///< Number of insertion attempts
        int acceptedInsertions{0};   ///< Number of accepted insertions
        int deletionAttempts{0};     ///< Number of deletion attempts
        int acceptedDeletions{0};    ///< Number of accepted deletions

        /**
         * @brief Get overall acceptance rate
         */
        double getAcceptanceRate() const {
            return (totalMoves > 0) ? static_cast<double>(acceptedMoves) / totalMoves : 0.0;
        }

        /**
         * @brief Get insertion acceptance rate
         */
        double getInsertionRate() const {
            return (insertionAttempts > 0) ? static_cast<double>(acceptedInsertions) / insertionAttempts : 0.0;
        }

        /**
         * @brief Get deletion acceptance rate
         */
        double getDeletionRate() const {
            return (deletionAttempts > 0) ? static_cast<double>(acceptedDeletions) / deletionAttempts : 0.0;
        }

        /**
         * @brief Reset all statistics
         */
        void reset() {
            totalMoves = 0;
            acceptedMoves = 0;
            insertionAttempts = 0;
            acceptedInsertions = 0;
            deletionAttempts = 0;
            acceptedDeletions = 0;
        }
    } stats;

    // Physical constants
    static constexpr float BOLTZMANN = 0.00831446f;  ///< Boltzmann constant [kJ/(mol·K)]
    static constexpr float MOLES_TO_MOLECULES = 0.0006023f;  
    static constexpr float MOLECULES_TO_MOLES = 1660.539f;   

    /**
     * @brief Set temperature and calculate beta
     */
    void setTemperature(float temperature) {
        if (temperature <= 0.0f) {
            throw std::invalid_argument("Temperature must be positive");
        }
        beta = 1.0f / (BOLTZMANN * temperature);
    }

    /**
     * @brief Get temperature from beta
     */
    float getTemperature() const {
        return (beta > 0.0f) ? 1.0f / (BOLTZMANN * beta) : 0.0f;
    }

    /**
     * @brief Set box dimensions
     */
    void setBox(float x, float y, float z) {
        box[0] = x;
        box[1] = y;
        box[2] = z;
        volume = x * y * z;
    }

    /**
     * @brief Get box volume
     */
    float getVolume() const {
        return volume;
    }
};

/**
 * @brief Force field parameters for Monte Carlo simulation
 */
struct MCForceField {
    int numTotalTypes;    ///< Total number of atom types in system
    int numMovementTypes; ///< Number of atom types in movement molecules

    /// @brief LJ sigma parameters [nm] - stored as matrix
    std::vector<float> ljSigma;   
    /// @brief LJ epsilon parameters [kJ/mol] - stored as matrix
    std::vector<float> ljEps;     

    /**
     * @brief Initialize force field for given number of types
     */
    void initialize(int total_types, int movement_types = 0) {
        numTotalTypes = total_types;
        numMovementTypes = movement_types;
        size_t matrix_size = total_types * total_types;
        ljSigma.resize(matrix_size, 0.0f);
        ljEps.resize(matrix_size, 0.0f);
    }

    /**
     * @brief Set LJ parameters for a type pair
     */
    void setLJParams(int type1, int type2, float sigma, float epsilon) {
        if (type1 >= numTotalTypes || type2 >= numTotalTypes) {
            throw std::out_of_range("Type index out of range");
        }
        int index = type1 * numTotalTypes + type2;
        ljSigma[index] = sigma;
        ljEps[index] = epsilon;
        
        // Set symmetric entry
        if (type1 != type2) {
            int sym_index = type2 * numTotalTypes + type1;
            ljSigma[sym_index] = sigma;
            ljEps[sym_index] = epsilon;
        }
    }

    /**
     * @brief Get LJ parameters for a type pair
     */
    std::pair<float, float> getLJParams(int type1, int type2) const {
        if (type1 >= numTotalTypes || type2 >= numTotalTypes) {
            throw std::out_of_range("Type index out of range");
        }
        int index = type1 * numTotalTypes + type2;
        return std::make_pair(ljSigma[index], ljEps[index]);
    }
};

/**
 * @brief Basic atomic properties for Monte Carlo simulation
 */
struct MCAtom {
    float x, y, z;      ///< Position (nm)
    float charge;       ///< Charge (e)
    int   type;        ///< Type index

    MCAtom() : x(0.0f), y(0.0f), z(0.0f), charge(0.0f), type(-1) {}
    MCAtom(float x_, float y_, float z_, float charge_, int type_) 
        : x(x_), y(y_), z(z_), charge(charge_), type(type_) {}

    /**
     * @brief Set position
     */
    void setPosition(float x_, float y_, float z_) {
        x = x_;
        y = y_;
        z = z_;
    }

    /**
     * @brief Get position as array
     */
    std::array<float, 3> getPosition() const {
        return {x, y, z};
    }
};

/**
 * @brief Molecular unit for GCMC simulation
 */
struct MCResidue {
    // Basic properties
    int   atomStart;    ///< Starting index in global atom array
    int   atomCount;    ///< Number of atoms in this residue
    bool  active;       ///< Whether this residue is currently in use
    bool  fixed;        ///< Whether this residue can be moved
    float center[3];    ///< Geometric center (nm)

    // Energy components
    float energy_vdw;   ///< Lennard-Jones energy (kJ/mol)
    float energy_elec;  ///< Coulomb energy (kJ/mol)
        
    // GCMC parameters
    float concentration; ///< Target concentration for GCMC (mol/L)
    float chemPot;      ///< Excess chemical potential (kJ/mol)
    int   type;         ///< Residue type index
    float radius;       ///< Approximate radius (nm)

    MCResidue() : atomStart(-1), atomCount(0), active(false), fixed(false),
                  center{0.0f, 0.0f, 0.0f}, energy_vdw(0.0f), energy_elec(0.0f),
                  concentration(0.0f), chemPot(0.0f), type(-1), radius(0.0f) {}

    /**
     * @brief Get total energy
     */
    float getTotalEnergy() const {
        return energy_vdw + energy_elec;
    }

    /**
     * @brief Set center position
     */
    void setCenter(float x, float y, float z) {
        center[0] = x;
        center[1] = y;
        center[2] = z;
    }

    /**
     * @brief Get center as array
     */
    std::array<float, 3> getCenter() const {
        return {center[0], center[1], center[2]};
    }

    /**
     * @brief Check if residue is valid
     */
    bool isValid() const {
        return atomStart >= 0 && atomCount > 0 && type >= 0;
    }
};

/**
 * @brief Information about movement residues in the system
 */
struct MCMovementResidueInfo {
    int startIndex;         ///< Starting index in global residue array
    int activeCount;        ///< Number of active movement residues
    int totalCount;         ///< Total number of movement residues
    std::string resName;    ///< Residue name for this movement group

    MCMovementResidueInfo() : startIndex(-1), activeCount(0), totalCount(0) {}
    MCMovementResidueInfo(int start, int active, int total, const std::string& name)
        : startIndex(start), activeCount(active), totalCount(total), resName(name) {}

    /**
     * @brief Get fraction of active residues
     */
    double getActiveFraction() const {
        return (totalCount > 0) ? static_cast<double>(activeCount) / totalCount : 0.0;
    }
};

/**
 * @brief Current state of the Monte Carlo system
 * This is the core state container that manages all system data
 */
class MCStateCore : public common::IValidatable {
public:
    MCStateCore() : activeAtomCount(0), activeResidueCount(0) {}

    // IValidatable interface
    bool is_valid() const override {
        // Check basic consistency
        if (activeAtomCount < 0 || activeResidueCount < 0) return false;
        if (activeAtomCount > static_cast<int>(atoms.size())) return false;
        if (activeResidueCount > static_cast<int>(residues.size())) return false;
        
        // Check residue-atom consistency
        for (int i = 0; i < activeResidueCount; ++i) {
            const auto& residue = residues[i];
            if (!residue.active) continue;
            if (residue.atomStart < 0 || residue.atomCount <= 0) return false;
            if (residue.atomStart + residue.atomCount > activeAtomCount) return false;
        }
        
        return true;
    }

    // Core data arrays
    std::vector<MCAtom>    atoms;      ///< Global atom array
    std::vector<MCResidue> residues;   ///< Global residue array

    // Type mappings
    TypeMaps residueTypes;  ///< Residue type mappings
    TypeMaps atomTypes;     ///< Atom type mappings

    // Movement molecule info
    std::vector<MCMovementResidueInfo> movementResidues;  ///< Info for movement residues
    std::vector<int> movementAtomTypes;  ///< Atom types from movement molecules
    int numMovementAtomTypes;  ///< Number of atom types from movement molecules

    // Active counts for swap-and-pop management
    int activeAtomCount;     ///< Current active atom count
    int activeResidueCount;  ///< Current active residue count

    // System info and force field
    MCInfo info;            ///< System information
    MCForceField forcefield;  ///< Force field parameters

    // Energy components
    struct EwaldEnergy {
        double real_space{0.0};     ///< Real space energy component (kJ/mol)
        double reciprocal{0.0};     ///< Reciprocal space energy component (kJ/mol)
        double self{0.0};           ///< Self-energy correction term (kJ/mol)
        double total{0.0};          ///< Total Ewald energy (kJ/mol)

        void reset() {
            real_space = 0.0;
            reciprocal = 0.0;
            self = 0.0;
            total = 0.0;
        }

        void updateTotal() {
            total = real_space + reciprocal + self;
        }
    } ewald_energy;

    // State management methods
    void reserve(int max_atoms, int max_residues) {
        atoms.reserve(max_atoms);
        residues.reserve(max_residues);
        info.maxAtoms = max_atoms;
        info.maxResidues = max_residues;
    }

    void clear() {
        atoms.clear();
        residues.clear();
        residueTypes.clear();
        atomTypes.clear();
        movementResidues.clear();
        movementAtomTypes.clear();
        activeAtomCount = 0;
        activeResidueCount = 0;
        numMovementAtomTypes = 0;
        ewald_energy.reset();
        info.stats.reset();
    }

    // Getters
    int getActiveAtomCount() const { return activeAtomCount; }
    int getActiveResidueCount() const { return activeResidueCount; }
    int getTotalAtomCount() const { return static_cast<int>(atoms.size()); }
    int getTotalResidueCount() const { return static_cast<int>(residues.size()); }

    // Energy calculation utilities
    double getTotalSystemEnergy() const {
        double total = ewald_energy.total;
        for (int i = 0; i < activeResidueCount; ++i) {
            if (residues[i].active) {
                total += residues[i].getTotalEnergy();
            }
        }
        return total;
    }
};

} // namespace montecarlo
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MONTECARLO_STATE_CORE_HPP 