#pragma once
#include "MCComposite.hpp"
#include "model/ModelModule.hpp"
#include <vector>
#include <string>

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief Movement state builder
 * 
 * Constructs the final Monte Carlo state by combining reindexed existing residues
 * with new movement molecules and their inactive copies.
 */
class MCMovementBuilder {
public:
    MCMovementBuilder() = default;
    ~MCMovementBuilder() = default;

    // Disable copy operations
    MCMovementBuilder(const MCMovementBuilder&) = delete;
    MCMovementBuilder& operator=(const MCMovementBuilder&) = delete;

    // Enable move operations
    MCMovementBuilder(MCMovementBuilder&&) = default;
    MCMovementBuilder& operator=(MCMovementBuilder&&) = default;

    /**
     * @brief Build final Monte Carlo state
     * 
     * Combines all components into the final state including:
     * - Unmatched residues (non-movement molecules)
     * - Matched residues (existing movement molecules)
     * - New inactive copies for movement molecules
     * 
     * @param newState State to populate
     * @param molecules Movement molecule information
     * @param processedResidueNames Processed residue names from type collector
     * @param reindexedResidues Reindexed existing residues
     * @param reindexedAtoms Reindexed atoms per residue
     */
    void buildFinalState(
        model::MCState& newState,
        const std::vector<MovementMolecularInfo>& molecules,
        const std::vector<std::string>& processedResidueNames,
        const std::vector<model::MCResidue>& reindexedResidues,
        const std::vector<std::vector<model::MCAtom>>& reindexedAtoms);

private:
    /**
     * @brief Match residues to movement molecules
     * 
     * @param reindexedResidues All reindexed residues
     * @param reindexedAtoms All reindexed atoms
     * @param processedResidueNames Target residue names to match
     * @param newState State containing type maps for name lookup
     * @return Tuple of (matching residues per molecule, matching atoms per molecule, other residues, other atoms)
     */
    std::tuple<
        std::vector<std::vector<model::MCResidue>>,
        std::vector<std::vector<std::vector<model::MCAtom>>>,
        std::vector<model::MCResidue>,
        std::vector<std::vector<model::MCAtom>>
    > matchResidues(
        const std::vector<model::MCResidue>& reindexedResidues,
        const std::vector<std::vector<model::MCAtom>>& reindexedAtoms,
        const std::vector<std::string>& processedResidueNames,
        const model::MCState& newState) const;

    /**
     * @brief Add unmatched residues to state
     * 
     * @param newState State to modify
     * @param otherResidues Unmatched residues
     * @param otherAtoms Atoms for unmatched residues
     * @param newAtomStart Current atom start position
     * @param newResIdx Current residue index
     */
    void addUnmatchedResidues(
        model::MCState& newState,
        const std::vector<model::MCResidue>& otherResidues,
        const std::vector<std::vector<model::MCAtom>>& otherAtoms,
        int& newAtomStart,
        int& newResIdx) const;

    /**
     * @brief Add movement molecule groups to state
     * 
     * @param newState State to modify
     * @param molecules Movement molecule information
     * @param processedResidueNames Processed residue names
     * @param matchingResidues Matched residues per molecule
     * @param matchingAtoms Matched atoms per molecule
     * @param newAtomStart Current atom start position
     * @param newResIdx Current residue index
     */
    void addMovementMoleculeGroups(
        model::MCState& newState,
        const std::vector<MovementMolecularInfo>& molecules,
        const std::vector<std::string>& processedResidueNames,
        const std::vector<std::vector<model::MCResidue>>& matchingResidues,
        const std::vector<std::vector<std::vector<model::MCAtom>>>& matchingAtoms,
        int& newAtomStart,
        int& newResIdx);

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