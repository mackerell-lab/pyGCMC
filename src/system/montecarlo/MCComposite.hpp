#pragma once
#include <memory>
#include <vector>
#include <set>
#include "model/ModelModule.hpp"
#include "MCCore.hpp"
#include "MCGeometry.hpp"
#include "MCSwitching.hpp"
#include "MCInitializer.hpp"

namespace pygcmc {
namespace system {
namespace montecarlo {

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
 * @brief Composite Monte Carlo system
 *
 * Combines all Monte Carlo subsystems into a unified interface.
 * Provides all functionality needed for GCMC simulations.
 */
class MCComposite {
public:
    MCComposite();
    ~MCComposite() = default;

    // Disable copy operations
    MCComposite(const MCComposite&) = delete;
    MCComposite& operator=(const MCComposite&) = delete;

    // Enable move operations
    MCComposite(MCComposite&&) = default;
    MCComposite& operator=(MCComposite&&) = default;

    // --------------------------------------------------------
    // System initialization and setup
    // --------------------------------------------------------

    /**
     * @brief Initialize the Monte Carlo system with given parameters
     * @param info System parameters (box size, temperature, cutoff, etc.)
     */
    void initialize(const model::MCInfo& info);

    /**
     * @brief Set force field parameters directly
     * @param ff Pre-configured force field parameters
     */
    void setForceField(const model::MCForceField& ff);

    /**
     * @brief Initialize force field from CHARMM parameters
     * @param ff CHARMM force field object
     * @throw std::runtime_error If required parameters are missing
     */
    void initializeForceField(const model::ForceField& ff);

    /**
     * @brief Initialize from molecular system (PDB/PSF/TOP)
     * @param molecular Molecular system
     * @throw std::runtime_error If molecular system is invalid
     */
    void initializeFromMolecular(const std::shared_ptr<model::Molecular>& molecular);

    /**
     * @brief Add initial system configuration
     * @param resVec Array of residues to add
     * @param resCount Number of residues
     * @param atomVec Array of atoms to add
     * @param atomCount Number of atoms
     * @throw std::runtime_error If system capacity is exceeded
     */
    void addInitialResidues(const model::MCResidue* resVec, int resCount,
                           const model::MCAtom* atomVec, int atomCount);

    /**
     * @brief Add movement molecules to the system
     * @param molecules List of movement molecules and their copy counts
     */
    void addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules);

    // --------------------------------------------------------
    // Core GCMC operations
    // --------------------------------------------------------

    /**
     * @brief Insert a new residue
     * @param res Residue to insert
     * @param atoms Atoms in the residue
     * @return Index of inserted residue, -1 if failed
     */
    int insertResidue(const model::MCResidue& res, const model::MCAtom* atoms);

    /**
     * @brief Remove a residue
     * @param resIdx Index of residue to remove
     * @return true if successful, false otherwise
     */
    bool removeResidue(int resIdx);

    /**
     * @brief Translate a residue
     * @param resIdx Residue index
     * @param dx X displacement [nm]
     * @param dy Y displacement [nm]
     * @param dz Z displacement [nm]
     */
    void translateResidue(int resIdx, float dx, float dy, float dz);

    // --------------------------------------------------------
    // Energy calculations
    // --------------------------------------------------------

    /**
     * @brief Calculate non-bonded energy between residues
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

    // --------------------------------------------------------
    // Switching function management
    // --------------------------------------------------------

    /**
     * @brief Set or disable CHARMM-style smooth switching function
     * @param enable Whether to enable the switching function
     * @param r_on Inner cutoff radius (nm)
     * @param r_off Outer cutoff radius (nm)
     */
    void setSwitchingFunction(bool enable, float r_on = 1.0f, float r_off = 1.2f);

    /**
     * @brief Calculate the switching function value at the given distance
     * @param r Distance (nm)
     * @return Switching function value [0,1]
     */
    float calculateSwitchingFunction(float r) const;

    // --------------------------------------------------------
    // State access
    // --------------------------------------------------------

    /**
     * @brief Get current system state
     * @return Const reference to current state
     */
    const model::MCState& getState() const { return state_; }

    /**
     * @brief Get modifiable system state
     * @return Modifiable reference to current state
     */
    model::MCState& getState() { return state_; }

    /**
     * @brief Get current type mappings
     * @return Current atom type mapping system
     */
    const model::TypeMaps& getTypeMaps() const { return state_.atomTypes; }

    /**
     * @brief Get number of active residues
     * @return Current active residue count
     */
    int getActiveResidueCount() const { return state_.activeResidueCount; }

    /**
     * @brief Get number of active atoms
     * @return Current active atom count
     */
    int getActiveAtomCount() const { return state_.activeAtomCount; }

    // --------------------------------------------------------
    // Parameter validation
    // --------------------------------------------------------

    /**
     * @brief Check if all atoms have valid topology and force field parameters
     * @param ff Force field to check against
     * @param molecular Molecular system to check
     * @throws std::runtime_error if any parameters are missing
     */
    void validateParameters(const model::ForceField& ff, const std::shared_ptr<model::Molecular>& molecular);

private:
    model::MCState state_;              ///< Current Monte Carlo state
    MCCore core_;                       ///< Core GCMC operations
    MCGeometry geometry_;               ///< Geometric calculations
    MCSwitching switching_;             ///< Switching function management
    MCInitializer initializer_;         ///< System initialization

    std::shared_ptr<model::Molecular> molecular_;  ///< Stored molecular system for validation
};

} // namespace montecarlo
} // namespace system
} // namespace pygcmc
