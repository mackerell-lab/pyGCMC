#pragma once
#include "model/montecarlo.hpp"

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief Geometric calculations for Monte Carlo system
 * 
 * Handles periodic boundary conditions and distance calculations
 * with proper minimum image convention.
 */
class MCGeometry {
public:
    MCGeometry() = default;
    ~MCGeometry() = default;

    // Disable copy operations
    MCGeometry(const MCGeometry&) = delete;
    MCGeometry& operator=(const MCGeometry&) = delete;

    // Enable move operations
    MCGeometry(MCGeometry&&) = default;
    MCGeometry& operator=(MCGeometry&&) = default;

    /**
     * @brief Calculate minimum image squared distance
     * 
     * Uses minimum image convention for periodic boundary conditions
     * consistent with CHARMM convention.
     * 
     * @param state Monte Carlo state containing box dimensions
     * @param dx,dy,dz Coordinate differences in nm
     * @return Squared distance in nm² for use in energy calculations
     */
    float getMinImageDistSqr(const model::MCState& state, float dx, float dy, float dz) const;

    /**
     * @brief Apply periodic boundary conditions to coordinates
     * 
     * @param state Monte Carlo state containing box dimensions
     * @param x,y,z Coordinates in nm (converted from PDB Å)
     * @note Uses box dimensions from PDB CRYST1 record (converted to nm)
     */
    void applyPBC(const model::MCState& state, float& x, float& y, float& z) const;

    /**
     * @brief Update residue geometric center
     * 
     * Calculates the geometric center of a residue based on atomic coordinates
     * 
     * @param state Monte Carlo state
     * @param res Residue to update
     */
    void updateGeometricCenter(const model::MCState& state, model::MCResidue& res) const;

    /**
     * @brief Calculate distance between two points with PBC
     * 
     * @param state Monte Carlo state containing box dimensions
     * @param x1,y1,z1 First point coordinates [nm]
     * @param x2,y2,z2 Second point coordinates [nm]
     * @return Distance in nm
     */
    float getDistance(const model::MCState& state, 
                     float x1, float y1, float z1,
                     float x2, float y2, float z2) const;

    /**
     * @brief Calculate squared distance between two points with PBC
     * 
     * @param state Monte Carlo state containing box dimensions
     * @param x1,y1,z1 First point coordinates [nm]
     * @param x2,y2,z2 Second point coordinates [nm]
     * @return Squared distance in nm²
     */
    float getDistanceSquared(const model::MCState& state, 
                            float x1, float y1, float z1,
                            float x2, float y2, float z2) const;
};

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 