#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_MERGER_HPP
#define PYGCMC_SYSTEM_MOLECULAR_MERGER_HPP

#include <memory>
#include <vector>
#include "model/molecular.hpp"
#include "model/topology.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

/**
 * @brief Merger for combining multiple topology files into a single molecular object
 */
class MolecularMerger {
public:
    MolecularMerger() = default;
    ~MolecularMerger() = default;

    /**
     * @brief Merge multiple topologies into a molecular object
     * @param molecular Target molecular object to merge into
     * @param topologies Vector of topology objects to merge
     */
    void mergeTopologies(
        std::shared_ptr<model::Molecular>& molecular,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);

private:
    /**
     * @brief Copy atoms from topology with offset adjustment
     * @param molecular Target molecular object
     * @param topology Source topology
     * @param atom_offset Offset to apply to atom IDs
     * @param residue_offset Offset to apply to residue IDs
     */
    void copyAtomsWithOffset(
        std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Topology>& topology,
        size_t atom_offset,
        size_t residue_offset);

    /**
     * @brief Copy residues from topology with offset adjustment
     * @param molecular Target molecular object
     * @param topology Source topology
     * @param atom_offset Offset to apply to atom IDs
     * @param residue_offset Offset to apply to residue IDs
     */
    void copyResiduesWithOffset(
        std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Topology>& topology,
        size_t atom_offset,
        size_t residue_offset);

    /**
     * @brief Copy bonding information with offset adjustment
     * @param molecular Target molecular object
     * @param topology Source topology
     * @param atom_offset Offset to apply to atom IDs
     */
    void copyBondingWithOffset(
        std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Topology>& topology,
        size_t atom_offset);

    /**
     * @brief Copy CMAP information with offset adjustment
     * @param molecular Target molecular object
     * @param topology Source topology
     * @param atom_offset Offset to apply to atom IDs
     * @param existing_cmaps Previously existing CMAPs to preserve
     */
    void copyCMAPWithOffset(
        std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Topology>& topology,
        size_t atom_offset,
        const std::vector<model::TopologyCmap>& existing_cmaps);
};

} // namespace molecular
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_MERGER_HPP 