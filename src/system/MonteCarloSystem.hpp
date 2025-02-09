/**
 * @file MonteCarloSystem.hpp
 * @brief Main class for managing Monte Carlo simulations
 * 
 * This class provides the core functionality for Grand Canonical Monte Carlo (GCMC)
 * simulations, including system initialization, molecule management, and energy calculations.
 */

#pragma once
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include "model/montecarlo.hpp"
#include "model/molecular.hpp"
#include "model/forcefield.hpp"

namespace pygcmc {
namespace system {

/**
 * @brief Main class for Monte Carlo simulation management
 * 
 * Handles all aspects of the Monte Carlo simulation, including:
 * - System initialization and setup
 * - Force field parameter management
 * - Molecule insertion, deletion, and movement
 * - Energy calculations
 * 
 * Uses swap-and-pop strategy for efficient memory management of active/inactive residues.
 */
class MonteCarloSystem {
public:
    // ------------------------------------------------------------
    // Constructor / Destructor
    // ------------------------------------------------------------
    MonteCarloSystem() = default;
    ~MonteCarloSystem() = default;

    // Disable copy operations to prevent accidental copies
    MonteCarloSystem(const MonteCarloSystem&) = delete;
    MonteCarloSystem& operator=(const MonteCarloSystem&) = delete;

    // Enable move operations for efficient container usage
    MonteCarloSystem(MonteCarloSystem&&) = default;
    MonteCarloSystem& operator=(MonteCarloSystem&&) = default;

    // ------------------------------------------------------------
    // System initialization and setup
    // ------------------------------------------------------------
    /**
     * @brief Initialize the Monte Carlo system with given parameters
     * 
     * @param info System parameters including box size, temperature, etc.
     * 
     * Allocates memory for atoms and residues based on maximum capacities
     * specified in the info structure.
     */
    void initialize(const model::MCInfo& info) {
        state.info = info;
        state.atoms.resize(info.maxAtoms);
        state.residues.resize(info.maxResidues);
    }

    /**
     * @brief Set force field parameters directly
     * 
     * @param ff Pre-configured force field parameters
     * 
     * Used when force field parameters have been pre-computed or loaded
     * from a previous simulation state.
     */
    void setForceField(const model::MCForceField& ff) {
        state.forcefield = ff;
    }

    /**
     * @brief Initialize force field from CHARMM parameters
     * 
     * @param ff CHARMM force field object
     * 
     * Processes CHARMM force field parameters to:
     * 1. Generate parameters for movement molecule interactions
     * 2. Apply combining rules
     * 3. Convert units to simulation units
     * 4. Handle NBFIX parameters
     * 
     * @throw std::runtime_error If parameters are missing
     */
    void initializeForceField(const model::ForceField& ff);

    /**
     * @brief Add initial system configuration
     * 
     * @param resVec Array of residues to add
     * @param resCount Number of residues
     * @param atomVec Array of atoms to add
     * @param atomCount Number of atoms
     * 
     * @throw std::runtime_error If system capacity is exceeded
     */
    void addInitialResidues(const model::MCResidue* resVec, int resCount,
                           const model::MCAtom* atomVec, int atomCount);

    /**
     * @brief Initialize from molecular system
     * 
     * @param molecular Molecular system to convert
     * 
     * Converts molecular system to Monte Carlo system:
     * 1. Converts coordinates and units
     * 2. Sets up type mappings
     * 3. Transfers molecular information
     * 
     * @throw std::runtime_error If molecular system is invalid
     */
    void initializeFromMolecular(const std::shared_ptr<model::Molecular>& molecular);

    /**
     * @brief Information for movement molecule initialization
     * 
     * Contains molecular structure and number of copies to pre-allocate
     * for each movement molecule type.
     */
    struct MovementMolecularInfo {
        std::shared_ptr<model::Molecular> molecular;  ///< Molecular structure
        int maxCopies;  ///< Number of copies to pre-allocate

        MovementMolecularInfo(std::shared_ptr<model::Molecular> mol, int max)
            : molecular(mol), maxCopies(max) {}
    };

    /**
     * @brief Add movement molecules to the system
     * 
     * @param molecules List of movement molecules and their copy counts
     * 
     * Sets up movement molecules for GCMC:
     * 1. Collects atom types
     * 2. Pre-allocates memory
     * 3. Initializes active and inactive copies
     * 
     * @note Order affects memory layout
     */
    void addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules);

    /**
     * @brief Get current type mappings
     * @return Current atom type mapping system
     */
    const model::TypeMaps& getTypeMaps() const { return state.atomTypes; }

    // ------------------------------------------------------------
    // Core GCMC operations
    // ------------------------------------------------------------
    /**
     * @brief Insert a new residue
     * 
     * @param res Residue to insert
     * @param atoms Atoms in the residue
     * @return Index of inserted residue, -1 if failed
     */
    int insertResidue(const model::MCResidue& res, const model::MCAtom* atoms);

    /**
     * @brief Remove a residue
     * 
     * @param resIdx Index of residue to remove
     * @return true if successful, false otherwise
     */
    bool removeResidue(int resIdx);

    /**
     * @brief Translate a residue
     * 
     * @param resIdx Residue index
     * @param dx X displacement [nm]
     * @param dy Y displacement [nm]
     * @param dz Z displacement [nm]
     */
    void translateResidue(int resIdx, float dx, float dy, float dz);

    /**
     * @brief Calculate non-bonded energy between residues
     * 
     * @param res1 First residue
     * @param res2 Second residue
     * @return Energy in kJ/mol
     */
    float calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const;

    /**
     * @brief Calculate total system energy
     * @return Total energy in kJ/mol
     */
    float calcTotalEnergy() const;

    // ------------------------------------------------------------
    // State access
    // ------------------------------------------------------------
    /**
     * @brief Get current system state
     * @return Const reference to current state
     */
    const model::MCState& getState() const { return state; }

    /**
     * @brief Get number of active residues
     * @return Current active residue count
     */
    int getActiveResidueCount() const { return state.activeResidueCount; }

    /**
     * @brief Get number of active atoms
     * @return Current active atom count
     */
    int getActiveAtomCount() const { return state.activeAtomCount; }

private:
    // ------------------------------------------------------------
    // Private helper functions
    // ------------------------------------------------------------
    /**
     * @brief Apply periodic boundary conditions
     * 
     * @param x X coordinate [nm]
     * @param y Y coordinate [nm]
     * @param z Z coordinate [nm]
     */
    void applyPBC(float& x, float& y, float& z) const;

    /**
     * @brief Calculate minimum image squared distance
     * 
     * @param dx X separation [nm]
     * @param dy Y separation [nm]
     * @param dz Z separation [nm]
     * @return Squared distance [nm²]
     */
    float getMinImageDistSqr(float dx, float dy, float dz) const;

    /**
     * @brief Update residue geometric center
     * @param res Residue to update
     */
    void updateGeometricCenter(model::MCResidue& res);

    // ------------------------------------------------------------
    // Member variables
    // ------------------------------------------------------------
    /// @brief Current system state
    model::MCState state;
};

} // namespace system
} // namespace pygcmc

