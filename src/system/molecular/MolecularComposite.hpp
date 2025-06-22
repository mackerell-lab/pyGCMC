#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_COMPOSITE_HPP
#define PYGCMC_SYSTEM_MOLECULAR_COMPOSITE_HPP

#include <memory>
#include <vector>
#include "model/ModelModule.hpp"
#include "MolecularCombiner.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

/**
 * @brief Composite interface for molecular system operations
 * 
 * This class provides a unified interface for all molecular system operations,
 * encapsulating the complexity of individual components.
 */
class MolecularComposite {
public:
    MolecularComposite();
    ~MolecularComposite() = default;

    /**
     * @brief Build molecular object from structure and topology
     * @param structure Structure data
     * @param topology Topology data
     * @return Combined molecular object
     */
    std::shared_ptr<model::Molecular> buildMolecular(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Build molecular object from structure and multiple topologies
     * @param structure Structure data
     * @param topologies Vector of topology data
     * @return Combined molecular object
     */
    std::shared_ptr<model::Molecular> buildMolecularMultiple(
        const std::shared_ptr<model::Structure>& structure,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);

    /**
     * @brief Get the current molecular object
     * @return Current molecular object
     */
    const std::shared_ptr<model::Molecular>& getCurrentMolecular() const { 
        return current_molecular_; 
    }

private:
    std::unique_ptr<MolecularCombiner> combiner_;
    std::shared_ptr<model::Molecular> current_molecular_;
};

} // namespace molecular
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_COMPOSITE_HPP 