#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_VALIDATOR_HPP
#define PYGCMC_SYSTEM_MOLECULAR_VALIDATOR_HPP

#include <memory>
#include <vector>
#include <string>
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

/**
 * @brief Validator for molecular system data consistency
 */
class MolecularValidator {
public:
    MolecularValidator() = default;
    ~MolecularValidator() = default;

    /**
     * @brief Verify atom types between structure and topology
     * @param pdb_res Residue from structure
     * @param top_res Residue from topology
     * @param topology Topology data
     */
    void verifyAtomTypes(
        const std::shared_ptr<model::Residue>& pdb_res,
        const model::TopologyResidue& top_res,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Validate overall combination consistency
     * @param structure Structure data
     * @param topology Topology data
     */
    void validateCombination(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief Validate multiple topology combination
     * @param structure Structure data
     * @param topologies Vector of topology data
     * @param total_atoms Expected total atoms
     * @param total_residues Expected total residues
     */
    void validateMultipleCombination(
        const std::shared_ptr<model::Structure>& structure,
        const std::vector<std::shared_ptr<model::Topology>>& topologies,
        size_t total_atoms,
        size_t total_residues);

private:
    /**
     * @brief Extract element from atom name or type
     * @param atom_name Atom name or type
     * @return Element string
     */
    std::string extractElement(const std::string& atom_name) const;

    /**
     * @brief Generate detailed error message for atom/residue mismatch
     * @param structure Structure data
     * @param topology Topology data
     * @param error_type Type of error
     * @return Detailed error message
     */
    std::string generateMismatchError(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology,
        const std::string& error_type) const;
};

} // namespace molecular
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_VALIDATOR_HPP 