#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_SPECIAL_HPP
#define PYGCMC_MODEL_TOPOLOGY_SPECIAL_HPP

#include "TopologyCore.hpp"
#include "TopologyValidation.hpp"
#include <stdexcept>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Special topology features management (groups, CMAP, exclusions)
 */
class TopologySpecialManager : public ValidationMixin {
public:
    TopologySpecialManager(TopologyStorage& storage) : ValidationMixin(storage), storage_(storage) {}

    /**
     * @brief Add a nonbonded exclusion between two atoms
     */
    void add_nonbonded_exclusion(int atom1, int atom2) {
        if (!are_valid_bond_atoms(atom1, atom2)) {
            throw std::invalid_argument("Invalid atom indices for exclusion");
        }
        
        storage_.exclusions[atom1].insert(atom2);
        storage_.exclusions[atom2].insert(atom1);
    }

    /**
     * @brief Add a group of atoms
     */
    void add_group(int id, const std::vector<int>& atoms, const std::string& type = "") {
        // Validate all atoms in the group
        for (int atom : atoms) {
            if (!is_valid_atom(atom)) {
                throw std::invalid_argument("Invalid atom index in group");
            }
        }
        
        TopologyGroup group;
        group.id = id;
        group.atoms = atoms;
        group.type = type;
        storage_.groups.push_back(group);
    }

    /**
     * @brief Add a CMAP term (8-atom CHARMM format)
     */
    void add_cmap(const std::array<int, 8>& atoms) {
        // Validate atoms (only first 5 are required to be valid for CHARMM format)
        for (int i = 0; i < 5; ++i) {
            if (!is_valid_atom(atoms[i])) {
                throw std::invalid_argument("Invalid atom index in CMAP");
            }
        }
        
        TopologyCmap cmap;
        cmap.atoms = atoms;
        cmap.function_type = 1;
        storage_.cmaps.push_back(cmap);
    }

    /**
     * @brief Add a CMAP term (5-atom GROMACS format)
     */
    void add_cmap(const std::array<int, 5>& atoms, int function_type = 1) {
        // Validate all atoms
        for (int atom : atoms) {
            if (!is_valid_atom(atom)) {
                throw std::invalid_argument("Invalid atom index in CMAP");
            }
        }
        
        // Convert 5-atom GROMACS format to 8-atom CHARMM format
        std::array<int, 8> charmm_atoms;
        for (int i = 0; i < 5; ++i) {
            charmm_atoms[i] = atoms[i];
        }
        for (int i = 5; i < 8; ++i) {
            charmm_atoms[i] = -1;
        }
        
        TopologyCmap cmap;
        cmap.atoms = charmm_atoms;
        cmap.function_type = function_type;
        storage_.cmaps.push_back(cmap);
    }

    // Getters
    const std::vector<TopologyGroup>& get_groups() const { return storage_.groups; }
    const std::vector<TopologyCmap>& get_cmaps() const { return storage_.cmaps; }
    const std::map<int, std::set<int>>& get_exclusions() const { return storage_.exclusions; }

    // Count methods
    size_t get_num_cmaps() const { return storage_.cmaps.size(); }
    size_t get_num_groups() const { return storage_.groups.size(); }

    // Check methods
    bool has_cmap() const {
        return !storage_.cmaps.empty();
    }

    bool has_group(int group_id) const {
        return std::any_of(storage_.groups.begin(), storage_.groups.end(),
            [group_id](const TopologyGroup& group) {
                return group.id == group_id;
            });
    }

    const TopologyGroup& get_group(int index) const {
        if (index < 0 || index >= static_cast<int>(storage_.groups.size())) {
            throw std::out_of_range("Invalid group index");
        }
        return storage_.groups[index];
    }

private:
    TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_SPECIAL_HPP 