#pragma once
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief Core Monte Carlo operations (insert, remove, translate)
 * 
 * Handles the fundamental GCMC operations with proper memory management
 * using swap-and-pop strategy for efficient residue handling.
 */
class MCCore {
public:
    MCCore() = default;
    ~MCCore() = default;

    // Disable copy operations
    MCCore(const MCCore&) = delete;
    MCCore& operator=(const MCCore&) = delete;

    // Enable move operations
    MCCore(MCCore&&) = default;
    MCCore& operator=(MCCore&&) = default;

    /**
     * @brief Insert a new residue into the system
     * 
     * @param state Monte Carlo state to operate on
     * @param res Residue to insert
     * @param atoms Atoms in the residue
     * @return Index of inserted residue, -1 if failed
     * 
     * @throw std::runtime_error If system capacity is exceeded
     */
    int insertResidue(model::MCState& state, const model::MCResidue& res, const model::MCAtom* atoms);

    /**
     * @brief Remove a residue from the system
     * 
     * Uses swap-and-pop strategy to maintain memory efficiency
     * 
     * @param state Monte Carlo state to operate on
     * @param resIdx Index of residue to remove
     * @return true if successful, false otherwise
     */
    bool removeResidue(model::MCState& state, int resIdx);

    /**
     * @brief Translate a residue by given displacement
     * 
     * Applies periodic boundary conditions automatically
     * 
     * @param state Monte Carlo state to operate on
     * @param resIdx Residue index
     * @param dx X displacement [nm]
     * @param dy Y displacement [nm] 
     * @param dz Z displacement [nm]
     */
    void translateResidue(model::MCState& state, int resIdx, float dx, float dy, float dz);

    /**
     * @brief Add initial residues to system
     * 
     * @param state Monte Carlo state to operate on
     * @param resVec Array of residues to add
     * @param resCount Number of residues
     * @param atomVec Array of atoms to add
     * @param atomCount Number of atoms
     * 
     * @throw std::runtime_error If initial system exceeds max capacity
     */
    void addInitialResidues(model::MCState& state, const model::MCResidue* resVec, int resCount,
                           const model::MCAtom* atomVec, int atomCount);

private:
    /**
     * @brief Apply periodic boundary conditions to coordinates
     * 
     * @param state Monte Carlo state containing box dimensions
     * @param x,y,z Coordinates to apply PBC [nm]
     */
    void applyPBC(const model::MCState& state, float& x, float& y, float& z) const;

    /**
     * @brief Update residue geometric center
     * 
     * @param state Monte Carlo state
     * @param res Residue to update
     */
    void updateGeometricCenter(const model::MCState& state, model::MCResidue& res) const;
};

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 