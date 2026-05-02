#pragma once
#include "MCComposite.hpp"
#include "model/ModelModule.hpp"
#include <vector>
#include <string>

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief Movement molecule type collector
 *
 * Collects and processes types from movement molecules for addMovementMolecules
 * operation. Handles type mapping expansion and validation.
 */
class MCMovementTypeCollector {
public:
    MCMovementTypeCollector() = default;
    ~MCMovementTypeCollector() = default;

    // Disable copy operations
    MCMovementTypeCollector(const MCMovementTypeCollector&) = delete;
    MCMovementTypeCollector& operator=(const MCMovementTypeCollector&) = delete;

    // Enable move operations
    MCMovementTypeCollector(MCMovementTypeCollector&&) = default;
    MCMovementTypeCollector& operator=(MCMovementTypeCollector&&) = default;

    /**
     * @brief Collect movement types from molecule list
     *
     * Processes the input molecules and extracts all unique residue types
     * and atom types, adding them to the provided type maps.
     *
     * @param molecules List of movement molecules
     * @param residueTypes Residue type map to expand
     * @param atomTypes Atom type map to expand
     * @return Vector of new movement atom type indices
     */
    std::vector<int> collectMovementTypes(
        const std::vector<MovementMolecularInfo>& molecules,
        model::TypeMaps& residueTypes,
        model::TypeMaps& atomTypes);

    /**
     * @brief Get processed residue names
     *
     * @return Vector of processed residue names (trimmed and uppercased)
     */
    const std::vector<std::string>& getProcessedResidueNames() const {
        return processedResidueNames_;
    }

private:
    std::vector<std::string> processedResidueNames_;  ///< Processed residue names for later use

    /**
     * @brief Process residue name (trim and uppercase)
     *
     * @param rawName Raw residue name from molecular data
     * @return Processed name (trimmed and uppercased)
     */
    std::string processResidueName(const std::string& rawName) const;

    /**
     * @brief Trim whitespace from string
     *
     * @param s Input string
     * @return Trimmed string
     */
    std::string trim(const std::string& s) const;
};

} // namespace montecarlo
} // namespace system
} // namespace pygcmc
