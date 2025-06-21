#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_COMBINER_HPP
#define PYGCMC_SYSTEM_MOLECULAR_COMBINER_HPP

#include <memory>
#include <vector>
#include "model/molecular.hpp"
#include "model/structure.hpp"
#include "model/topology.hpp"
#include "MolecularValidator.hpp"
#include "MolecularMatcher.hpp"
#include "MolecularMerger.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

/**
 * @brief Combiner for structure and topology data into molecular objects
 */
class MolecularCombiner {
public:
    MolecularCombiner();
    ~MolecularCombiner() = default;

    /**
     * @brief Combine single structure and topology into molecular object
     * @param structure Structure data
     * @param topology Topology data
     * @return Combined molecular object
     */
    std::shared_ptr<model::Molecular> combine(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Combine structure with multiple topologies
     * @param structure Structure data
     * @param topologies Vector of topology data
     * @return Combined molecular object
     */
    std::shared_ptr<model::Molecular> combineMultiple(
        const std::shared_ptr<model::Structure>& structure,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);

private:
    std::unique_ptr<MolecularValidator> validator_;
    std::unique_ptr<MolecularMatcher> matcher_;
    std::unique_ptr<MolecularMerger> merger_;

    /**
     * @brief Copy basic structure data to molecular object
     * @param molecular Target molecular object
     * @param structure Source structure data
     */
    void copyStructureData(
        std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Structure>& structure);

    /**
     * @brief Copy topology data to molecular object
     * @param molecular Target molecular object
     * @param topology Source topology data
     */
    void copyTopologyData(
        std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Build lookup mappings for molecular object
     * @param molecular Target molecular object
     * @param topology Source topology data
     */
    void buildLookupMappings(
        std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Match residues between structure and topologies
     * @param structure Structure data
     * @param topologies Vector of topology data
     * @return Vector of matched topologies
     */
    std::vector<std::shared_ptr<model::Topology>> matchResidues(
        const std::shared_ptr<model::Structure>& structure,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);
};

} // namespace molecular
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_COMBINER_HPP 