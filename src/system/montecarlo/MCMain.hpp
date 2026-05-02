#pragma once
#include "MCComposite.hpp"
#include "model/ModelModule.hpp"
#include <memory>

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief Main class for Monte Carlo simulation management (Compatibility Layer)
 *
 * This is a compatibility layer that maintains the original MonteCarloSystem API
 * while internally using the new modular architecture. All original functionality
 * is preserved while gaining the benefits of the refactored design.
 */
class MCMain {
public:
    // ------------------------------------------------------------
    // Constructor / Destructor
    // ------------------------------------------------------------
    MCMain() = default;
    ~MCMain() = default;

    // Disable copy operations to prevent accidental copies
    MCMain(const MCMain&) = delete;
    MCMain& operator=(const MCMain&) = delete;

    // Enable move operations for efficient container usage
    MCMain(MCMain&&) = default;
    MCMain& operator=(MCMain&&) = default;

    // ------------------------------------------------------------
    // System initialization and setup
    // ------------------------------------------------------------

    /**
     * @brief Initialize the Monte Carlo system with given parameters
     * @param info System parameters including box size, temperature, cutoff distance
     */
    void initialize(const model::MCInfo& info) {
        impl_.initialize(info);
    }

    /**
     * @brief Set force field parameters directly
     * @param ff Pre-configured force field parameters
     */
    void setForceField(const model::MCForceField& ff) {
        impl_.setForceField(ff);
    }

    /**
     * @brief Initialize force field from CHARMM parameters
     * @param ff CHARMM force field object
     * @throw std::runtime_error If required parameters are missing
     */
    void initializeForceField(const model::ForceField& ff) {
        impl_.initializeForceField(ff);
    }

    /**
     * @brief Add initial system configuration
     * @param resVec Array of residues to add
     * @param resCount Number of residues
     * @param atomVec Array of atoms to add
     * @param atomCount Number of atoms
     * @throw std::runtime_error If system capacity is exceeded
     */
    void addInitialResidues(const model::MCResidue* resVec, int resCount,
                           const model::MCAtom* atomVec, int atomCount) {
        impl_.addInitialResidues(resVec, resCount, atomVec, atomCount);
    }

    /**
     * @brief Initialize from molecular system (PDB/PSF/TOP)
     * @param molecular Molecular system
     * @throw std::runtime_error If molecular system is invalid
     */
    void initializeFromMolecular(const std::shared_ptr<model::Molecular>& molecular) {
        impl_.initializeFromMolecular(molecular);
    }

    /**
     * @brief Add movement molecules to the system
     * @param molecules List of movement molecules and their copy counts
     */
    void addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules) {
        // Since we're in the same namespace, MovementMolecularInfo is the same type
        impl_.addMovementMolecules(molecules);
    }

    /**
     * @brief Get current type mappings
     * @return Current atom type mapping system
     */
    const model::TypeMaps& getTypeMaps() const {
        return impl_.getTypeMaps();
    }

    // ------------------------------------------------------------
    // Core GCMC operations
    // ------------------------------------------------------------

    /**
     * @brief Insert a new residue
     * @param res Residue to insert
     * @param atoms Atoms in the residue
     * @return Index of inserted residue, -1 if failed
     */
    int insertResidue(const model::MCResidue& res, const model::MCAtom* atoms) {
        return impl_.insertResidue(res, atoms);
    }

    /**
     * @brief Remove a residue
     * @param resIdx Index of residue to remove
     * @return true if successful, false otherwise
     */
    bool removeResidue(int resIdx) {
        return impl_.removeResidue(resIdx);
    }

    /**
     * @brief Translate a residue
     * @param resIdx Residue index
     * @param dx X displacement [nm]
     * @param dy Y displacement [nm]
     * @param dz Z displacement [nm]
     */
    void translateResidue(int resIdx, float dx, float dy, float dz) {
        impl_.translateResidue(resIdx, dx, dy, dz);
    }

    /**
     * @brief Calculate non-bonded energy between residues using CHARMM parameters
     * @param res1 First residue
     * @param res2 Second residue
     * @return Energy in kJ/mol (converted from CHARMM kcal/mol)
     */
    float calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const {
        return impl_.calcNonBondedEnergy(res1, res2);
    }

    /**
     * @brief Calculate total system energy using CHARMM parameters
     * @return Total energy in kJ/mol (converted from CHARMM kcal/mol)
     */
    float calcTotalEnergy() const {
        return impl_.calcTotalEnergy();
    }

    // ------------------------------------------------------------
    // State access
    // ------------------------------------------------------------

    /**
     * @brief Get current system state
     * @return Const reference to current state
     */
    const model::MCState& getState() const {
        return impl_.getState();
    }

    /**
     * @brief Get modifiable system state
     * @return Modifiable reference to current state
     */
    model::MCState& getState() {
        return impl_.getState();
    }

    /**
     * @brief Get number of active residues
     * @return Current active residue count
     */
    int getActiveResidueCount() const {
        return impl_.getActiveResidueCount();
    }

    /**
     * @brief Get number of active atoms
     * @return Current active atom count
     */
    int getActiveAtomCount() const {
        return impl_.getActiveAtomCount();
    }

    /**
     * @brief Check if all atoms have valid topology and force field parameters
     * @param ff Force field to check against
     * @param molecular Molecular system to check
     * @throws std::runtime_error if any parameters are missing
     */
    void validateParameters(const model::ForceField& ff, const std::shared_ptr<model::Molecular>& molecular) {
        impl_.validateParameters(ff, molecular);
    }

    /**
     * @brief Set or disable CHARMM-style smooth switching function
     * @param enable Whether to enable the switching function
     * @param r_on Inner cutoff radius (nm)
     * @param r_off Outer cutoff radius (nm)
     */
    void setSwitchingFunction(bool enable, float r_on = 1.0f, float r_off = 1.2f) {
        impl_.setSwitchingFunction(enable, r_on, r_off);
    }

    /**
     * @brief Calculate the switching function value at the given distance
     * @param r Distance (nm)
     * @return Switching function value [0,1]
     */
    float calculateSwitchingFunction(float r) const {
        return impl_.calculateSwitchingFunction(r);
    }

    /**
     * @brief Get whether the switching function is currently enabled
     * @return Whether the switching function is enabled
     */
    bool isUsingSwitchingFunction() const {
        return impl_.getState().info.use_switching;
    }

    /**
     * @brief Get inner cutoff radius
     * @return Inner cutoff radius (nm)
     */
    float getSwitchingROn() const {
        return impl_.getState().info.r_on;
    }

    /**
     * @brief Get outer cutoff radius
     * @return Outer cutoff radius (nm)
     */
    float getSwitchingROff() const {
        return impl_.getState().info.r_off;
    }

    /**
     * @brief Apply current switching function settings to external state object
     * @param externalState External state object to apply settings to
     */
    void applySwitchingToState(model::MCState& externalState) const {
        const auto& state = impl_.getState();
        externalState.info.use_switching = state.info.use_switching;
        externalState.info.r_on = state.info.r_on;
        externalState.info.r_off = state.info.r_off;
    }

private:
    MCComposite impl_;  ///< Internal implementation using new modular architecture
};

} // namespace montecarlo

// Export to parent namespace for backward compatibility
using MonteCarloSystem = montecarlo::MCMain;
using MovementMolecularInfo = montecarlo::MovementMolecularInfo;

} // namespace system
} // namespace pygcmc
