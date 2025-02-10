/**
 * @file MonteCarloSystem.hpp
 * @brief Main class for managing Monte Carlo simulations
 * 
 * This class provides the core functionality for Grand Canonical Monte Carlo (GCMC)
 * simulations, including system initialization, molecule management, and energy calculations.
 * 
 * Handles unit conversions between PDB/CHARMM and internal units:
 * - Coordinates: PDB Å -> internal nm
 * - Energies: CHARMM kcal/mol -> internal kJ/mol
 * - Charges: CHARMM partial charges (e) used directly
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
 * - System initialization and setup from PDB/CHARMM inputs
 * - Force field parameter conversion and management
 * - Molecule insertion, deletion, and movement with PBC
 * - Energy calculations using CHARMM parameters
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
     * @param info System parameters including:
     *        - Box size (nm, converted from PDB CRYST1 record)
     *        - Temperature (K)
     *        - Cutoff distance (nm)
     * 
     * Allocates memory for atoms and residues based on maximum capacities.
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
     * @param ff CHARMM force field object containing:
     *        - Atom types and charges
     *        - LJ parameters (Rmin/2 in Å, epsilon in kcal/mol)
     *        - NBFIX parameters if available
     * 
     * Processes CHARMM parameters:
     * 1. Converts units (Å -> nm, kcal/mol -> kJ/mol)
     * 2. Transforms Rmin/2 to sigma
     * 3. Applies combining rules or NBFIX
     * 4. Organizes parameters for efficient access
     * 
     * @throw std::runtime_error If required parameters are missing
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
     * @brief Initialize from molecular system (PDB/PSF/TOP)
     * 
     * @param molecular Molecular system containing:
     *        - Atomic coordinates (Å)
     *        - CHARMM atom types and charges
     *        - Residue definitions
     *        - Box dimensions (Å)
     * 
     * Performs:
     * 1. Coordinate conversion (Å -> nm)
     * 2. Box dimension conversion
     * 3. Type mapping setup
     * 4. Memory allocation and data transfer
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
     * @brief Calculate non-bonded energy between residues using CHARMM parameters
     * 
     * Includes:
     * 1. Lennard-Jones with converted CHARMM parameters:
     *    - Uses sigma-epsilon form converted from CHARMM Rmin/2-epsilon
     *    - Energy in kJ/mol (converted from CHARMM kcal/mol)
     * 2. Coulomb with CHARMM partial charges
     * 3. Periodic boundary conditions
     * 4. Standard CHARMM cutoff scheme
     * 
     * @param res1 First residue
     * @param res2 Second residue
     * @return Energy in kJ/mol (converted from CHARMM kcal/mol)
     */
    float calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const;

    /**
     * @brief Calculate total system energy using CHARMM parameters
     * 
     * Computes:
     * 1. All non-bonded interactions using converted CHARMM parameters
     * 2. Applies periodic boundary conditions
     * 3. Uses cutoff-based neighbor lists
     * 
     * @return Total energy in kJ/mol (converted from CHARMM kcal/mol)
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

    /**
     * @brief Check if all atoms have valid topology and force field parameters
     * @param ff Force field to check against
     * @param molecular Molecular system to check
     * @throws std::runtime_error if any parameters are missing
     */
    void validateParameters(const model::ForceField& ff, const std::shared_ptr<model::Molecular>& molecular);

private:
    // ------------------------------------------------------------
    // Private helper functions
    // ------------------------------------------------------------
    /**
     * @brief Apply periodic boundary conditions to coordinates
     * 
     * @param x,y,z Coordinates in nm (converted from PDB Å)
     * @note Uses box dimensions from PDB CRYST1 record (converted to nm)
     */
    void applyPBC(float& x, float& y, float& z) const;

    /**
     * @brief Calculate minimum image squared distance
     * 
     * @param dx,dy,dz Coordinate differences in nm
     * @return Squared distance in nm² for use in energy calculations
     * @note Consistent with CHARMM's minimum image convention
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
    model::MCState state;  ///< Current state of the Monte Carlo system
    std::shared_ptr<model::Molecular> molecular;  ///< Stored molecular system for parameter validation
};

} // namespace system
} // namespace pygcmc

