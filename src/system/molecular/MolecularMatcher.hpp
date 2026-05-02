#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_MATCHER_HPP
#define PYGCMC_SYSTEM_MOLECULAR_MATCHER_HPP

#include <memory>
#include <vector>
#include <string>
#include <set>
#include <map>
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

// Forward declaration
class MolecularValidator;

/**
 * @brief Matcher for residue sequences between structure and topology
 */
class MolecularMatcher {
public:
    MolecularMatcher();
    ~MolecularMatcher() = default;

    /**
     * @brief Match residue sequence starting from a given index
     * @param pdb_residues Residues from structure
     * @param topology Topology data
     * @param start_idx Starting index for matching
     * @param matched_count Number of matched residues (output)
     * @return True if sequence matches successfully
     */
    bool matchResidueSequence(
        const std::vector<std::shared_ptr<model::Residue>>& pdb_residues,
        const std::shared_ptr<model::Topology>& topology,
        size_t start_idx,
        size_t& matched_count);

    /**
     * @brief Perform detailed residue matching with protein/non-protein classification
     * @param molecular Molecular object to validate against
     * @param topology Topology data
     */
    void performDetailedMatching(
        const std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Topology>& topology);

private:
    std::unique_ptr<MolecularValidator> validator_;

    // Standard amino acid list
    static const std::set<std::string> standard_amino_acids_;

    /**
     * @brief Build residue connection relationships from bonds
     * @param topology Topology data
     * @return Map of residue connections
     */
    std::map<int, std::set<int>> buildResidueConnections(
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Find connected residue chains using BFS
     * @param residue_connections Residue connection map
     * @param num_residues Total number of residues
     * @return Vector of chains (each chain is a set of residue IDs)
     */
    std::vector<std::set<int>> findConnectedChains(
        const std::map<int, std::set<int>>& residue_connections,
        size_t num_residues);

    /**
     * @brief Calculate amino acid count for each chain
     * @param chains Connected chains
     * @param topology Topology data
     * @return Map of chain to amino acid count
     */
    std::map<std::set<int>, int> calculateChainAACount(
        const std::vector<std::set<int>>& chains,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Check if a residue is a standard amino acid
     * @param resname Residue name
     * @return True if it's a standard amino acid
     */
    bool isStandardAminoAcid(const std::string& resname) const;
};

} // namespace molecular
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_MATCHER_HPP
