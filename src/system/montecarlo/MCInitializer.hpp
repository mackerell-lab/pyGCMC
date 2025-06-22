#pragma once
#include "model/ModelModule.hpp"
#include <vector>
#include <memory>

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief Monte Carlo system initializer
 * 
 * Handles initialization of Monte Carlo system from molecular data,
 * including unit conversions and parameter validation.
 */
class MCInitializer {
public:
    MCInitializer() = default;
    ~MCInitializer() = default;

    // Disable copy operations
    MCInitializer(const MCInitializer&) = delete;
    MCInitializer& operator=(const MCInitializer&) = delete;

    // Enable move operations
    MCInitializer(MCInitializer&&) = default;
    MCInitializer& operator=(MCInitializer&&) = default;

    /**
     * @brief Initialize Monte Carlo state from molecular system (PDB/PSF/TOP)
     * 
     * @param state Monte Carlo state to initialize
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
     * @throw std::runtime_error If molecular system is invalid or exceeds capacity
     */
    void initializeFromMolecular(model::MCState& state, const std::shared_ptr<model::Molecular>& molecular);

    /**
     * @brief Initialize force field from CHARMM parameters
     * 
     * @param state Monte Carlo state to modify
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
    void initializeForceField(model::MCState& state, const model::ForceField& ff);

    /**
     * @brief Validate parameters consistency
     * 
     * @param ff Force field to check against
     * @param molecular Molecular system to check
     * @throw std::runtime_error if any parameters are missing
     */
    void validateParameters(const model::ForceField& ff, const std::shared_ptr<model::Molecular>& molecular);

private:
    // Unit conversion constants
    static constexpr float ANGSTROM_TO_NM = 0.1f;    ///< 1 Å = 0.1 nm
    static constexpr float KCAL_TO_KJ = 4.184f;      ///< 1 kcal/mol = 4.184 kJ/mol

    /**
     * @brief Convert molecular residues to Monte Carlo residues
     * 
     * @param state Monte Carlo state to populate
     * @param molecular Source molecular system
     * @return Pair of (converted residues, converted atoms)
     */
    std::pair<std::vector<model::MCResidue>, std::vector<model::MCAtom>>
    convertMolecularData(model::MCState& state, const std::shared_ptr<model::Molecular>& molecular);
};

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 