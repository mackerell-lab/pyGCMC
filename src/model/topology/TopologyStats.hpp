#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_STATS_HPP
#define PYGCMC_MODEL_TOPOLOGY_STATS_HPP

#include "TopologyCore.hpp"
#include <sstream>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Topology statistics structure
 */
struct TopologyStats {
    size_t num_atoms = 0;
    size_t num_residues = 0;
    size_t num_segments = 0;
    size_t num_bonds = 0;
    size_t num_angles = 0;
    size_t num_dihedrals = 0;
    size_t num_impropers = 0;
    size_t num_donors = 0;
    size_t num_acceptors = 0;
    size_t num_groups = 0;
    size_t num_cmaps = 0;
    size_t num_exclusions = 0;
};

/**
 * @brief Topology statistics generator
 */
class TopologyStatsGenerator {
public:
    TopologyStatsGenerator(const TopologyStorage& storage) : storage_(storage) {}

    /**
     * @brief Get comprehensive topology statistics
     */
    TopologyStats get_statistics() const {
        TopologyStats stats;
        stats.num_atoms = storage_.atoms.size();
        stats.num_residues = storage_.residues.size();
        stats.num_segments = storage_.segments.size();
        stats.num_bonds = storage_.bonds.size();
        stats.num_angles = storage_.angles.size();
        
        stats.num_dihedrals = 0;
        stats.num_impropers = 0;
        for (const auto& dih : storage_.dihedrals) {
            if (dih.improper) {
                stats.num_impropers++;
            } else {
                stats.num_dihedrals++;
            }
        }
        
        stats.num_donors = storage_.donors.size();
        stats.num_acceptors = storage_.acceptors.size();
        stats.num_groups = storage_.groups.size();
        stats.num_cmaps = storage_.cmaps.size();
        
        stats.num_exclusions = 0;
        for (const auto& pair : storage_.exclusions) {
            stats.num_exclusions += pair.second.size();
        }
        stats.num_exclusions /= 2; // Each exclusion is counted twice
        
        return stats;
    }

    /**
     * @brief Get human-readable topology summary
     */
    std::string get_summary() const {
        auto stats = get_statistics();
        std::stringstream ss;
        ss << "Topology Summary:\n";
        ss << "  Atoms: " << stats.num_atoms << "\n";
        ss << "  Residues: " << stats.num_residues << "\n";
        ss << "  Segments: " << stats.num_segments << "\n";
        ss << "  Bonds: " << stats.num_bonds << "\n";
        ss << "  Angles: " << stats.num_angles << "\n";
        ss << "  Dihedrals: " << stats.num_dihedrals << "\n";
        ss << "  Impropers: " << stats.num_impropers << "\n";
        if (stats.num_donors > 0) ss << "  H-bond donors: " << stats.num_donors << "\n";
        if (stats.num_acceptors > 0) ss << "  H-bond acceptors: " << stats.num_acceptors << "\n";
        if (stats.num_groups > 0) ss << "  Groups: " << stats.num_groups << "\n";
        if (stats.num_cmaps > 0) ss << "  CMAP terms: " << stats.num_cmaps << "\n";
        if (stats.num_exclusions > 0) ss << "  Exclusions: " << stats.num_exclusions << "\n";
        return ss.str();
    }

private:
    const TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_STATS_HPP 