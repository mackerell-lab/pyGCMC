#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_ANALYSIS_HPP
#define PYGCMC_MODEL_TOPOLOGY_ANALYSIS_HPP

#include "TopologyCore.hpp"
#include "TopologySearch.hpp"
#include <map>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Topology analysis utilities
 */
class TopologyAnalyzer {
public:
    TopologyAnalyzer(const TopologyStorage& storage) 
        : storage_(storage), searcher_(storage) {}

    /**
     * @brief Count atoms by type
     */
    std::map<std::string, int> count_atom_types() const {
        std::map<std::string, int> counts;
        for (const auto& atom : storage_.atoms) {
            counts[atom.type]++;
        }
        return counts;
    }

    /**
     * @brief Count residues by name
     */
    std::map<std::string, int> count_residue_types() const {
        std::map<std::string, int> counts;
        for (const auto& residue : storage_.residues) {
            counts[residue.name]++;
        }
        return counts;
    }

    /**
     * @brief Find shortest path between atoms (bond distance)
     */
    int get_bond_distance(int atom1, int atom2) const {
        if (atom1 == atom2) return 0;
        
        const size_t num_atoms = storage_.atoms.size();
        std::vector<bool> visited(num_atoms, false);
        std::vector<int> queue, distance(num_atoms, -1);
        
        queue.push_back(atom1);
        visited[atom1] = true;
        distance[atom1] = 0;
        
        for (size_t front = 0; front < queue.size(); ++front) {
            int current = queue[front];
            auto neighbors = searcher_.get_bonded_atoms(current);
            
            for (int neighbor : neighbors) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    distance[neighbor] = distance[current] + 1;
                    queue.push_back(neighbor);
                    
                    if (neighbor == atom2) {
                        return distance[neighbor];
                    }
                }
            }
        }
        
        return -1; // Not connected
    }

    /**
     * @brief Get connected components analysis
     */
    std::vector<std::vector<int>> get_connected_components() const {
        const size_t num_atoms = storage_.atoms.size();
        std::vector<bool> visited(num_atoms, false);
        std::vector<std::vector<int>> components;
        
        for (size_t i = 0; i < num_atoms; ++i) {
            if (!visited[i]) {
                std::vector<int> component;
                std::vector<int> stack = {static_cast<int>(i)};
                
                while (!stack.empty()) {
                    int current = stack.back();
                    stack.pop_back();
                    
                    if (!visited[current]) {
                        visited[current] = true;
                        component.push_back(current);
                        
                        auto neighbors = searcher_.get_bonded_atoms(current);
                        for (int neighbor : neighbors) {
                            if (!visited[neighbor]) {
                                stack.push_back(neighbor);
                            }
                        }
                    }
                }
                
                if (!component.empty()) {
                    components.push_back(std::move(component));
                }
            }
        }
        
        return components;
    }

private:
    const TopologyStorage& storage_;
    TopologySearcher searcher_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_ANALYSIS_HPP 