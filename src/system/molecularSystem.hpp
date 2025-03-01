#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP
#define PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP

#include <memory>
#include <vector>
#include "model/molecular.hpp"
#include "model/structure.hpp"
#include "model/topology.hpp"

namespace pygcmc {
namespace system {

/**
 * @brief Molecular system class for managing and operating on molecular data
 */
class MolecularSystem {
public:
    MolecularSystem() = default;
    ~MolecularSystem() = default;

    /**
     * @brief Combine Structure and Topology data into a Molecular object
     * @param structure Structure data
     * @param topology Topology data
     * @return Combined Molecular object
     */
    std::shared_ptr<model::Molecular> combine(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Combine multiple Structure and Topology data into a Molecular object
     * @param structure Structure data
     * @param topologies Topology data
     * @return Combined Molecular object
     */
    std::shared_ptr<model::Molecular> combine_multiple(
        const std::shared_ptr<model::Structure>& structure,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);

    /**
     * @brief Get the current Molecular object
     * @return Current Molecular object
     */
    const std::shared_ptr<model::Molecular>& get_molecular() const { return molecular_; }

private:
    std::shared_ptr<model::Molecular> molecular_;

    /**
     * @brief Check if residue sequence matches
     * @param pdb_residues Residue sequence in the structure
     * @param topology Topology data
     * @param start_idx Starting index
     * @param matched_count Number of matched residues
     * @return Whether it matches
     */
    bool match_residue_sequence(
        const std::vector<std::shared_ptr<model::Residue>>& pdb_residues,
        const std::shared_ptr<model::Topology>& topology,
        size_t start_idx,
        size_t& matched_count);

    /**
     * @brief Merge multiple topology files
     * @param molecular Molecular object
     * @param topologies Topology data
     */
    void merge_topologies(
        std::shared_ptr<model::Molecular>& molecular,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);

    /**
     * @brief Verify atom type matching
     * @param pdb_res Residue in the structure
     * @param top_res Residue in the topology
     * @param topology Topology data
     */
    void verify_atom_types(
        const std::shared_ptr<model::Residue>& pdb_res,
        const model::TopologyResidue& top_res,
        const std::shared_ptr<model::Topology>& topology);
};

} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP
