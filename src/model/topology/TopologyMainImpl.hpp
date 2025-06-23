#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_MAIN_IMPL_HPP
#define PYGCMC_MODEL_TOPOLOGY_MAIN_IMPL_HPP

#include "TopologyMain.hpp"
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

// Implementation of the only method that is not inlined in TopologyMain.hpp
inline bool Topology::has_cmap(const std::vector<int>& atoms) const {
    if (atoms.size() < 5) return false;
    return std::any_of(special_manager_.get_cmaps().begin(), special_manager_.get_cmaps().end(),
        [&atoms](const TopologyCmap& cmap) {
            for (size_t i = 0; i < 5 && i < atoms.size(); ++i) {
                if (cmap.atoms[i] != atoms[i]) return false;
            }
            return true;
        });
}

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_MAIN_IMPL_HPP 