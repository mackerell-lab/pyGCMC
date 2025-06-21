#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_MAIN_HPP
#define PYGCMC_SYSTEM_MOLECULAR_MAIN_HPP

#include <memory>
#include <vector>
#include "model/molecular.hpp"
#include "model/structure.hpp"
#include "model/topology.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

// Forward declaration
class MolecularComposite;

/**
 * @brief Main interface for molecular system operations
 * 
 * This class provides backward compatibility with the original MolecularSystem class
 * while using the new modular architecture internally.
 */
class MolecularMain {
public:
    MolecularMain();
    ~MolecularMain() = default;

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
    const std::shared_ptr<model::Molecular>& get_molecular() const;

private:
    std::unique_ptr<MolecularComposite> impl_;
};

} // namespace molecular
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_MAIN_HPP 