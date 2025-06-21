#pragma once
#include "model/montecarlo.hpp"
#include <vector>
#include <string>

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief Movement residue reindexer
 * 
 * Handles reindexing of existing residues and atoms when adding movement molecules.
 * This ensures type consistency across the expanded type maps.
 */
class MCMovementReindexer {
public:
    MCMovementReindexer() = default;
    ~MCMovementReindexer() = default;

    // Disable copy operations
    MCMovementReindexer(const MCMovementReindexer&) = delete;
    MCMovementReindexer& operator=(const MCMovementReindexer&) = delete;

    // Enable move operations
    MCMovementReindexer(MCMovementReindexer&&) = default;
    MCMovementReindexer& operator=(MCMovementReindexer&&) = default;

    /**
     * @brief Reindex existing residues and atoms
     * 
     * Takes the current system state and reindexes all residues and atoms
     * to use the new expanded type maps.
     * 
     * @param oldState Original state with old type maps
     * @param newResidueTypes New expanded residue type map
     * @param newAtomTypes New expanded atom type map
     * @return Pair of (reindexed residues, reindexed atoms per residue)
     */
    std::pair<std::vector<model::MCResidue>, std::vector<std::vector<model::MCAtom>>>
    reindexExistingResidues(const model::MCState& oldState,
                           model::TypeMaps& newResidueTypes,
                           model::TypeMaps& newAtomTypes);

private:
    /**
     * @brief Trim whitespace from string
     * 
     * @param s Input string
     * @return Trimmed string
     */
    std::string trim(const std::string& s) const;

    /**
     * @brief Process residue name (trim and uppercase)
     * 
     * @param rawName Raw residue name
     * @return Processed name (trimmed and uppercased)
     */
    std::string processResidueName(const std::string& rawName) const;
};

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 