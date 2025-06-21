#include "MolecularCombiner.hpp"
#include "MolecularValidator.hpp"
#include "MolecularMatcher.hpp"
#include "MolecularMerger.hpp"
#include <stdexcept>
#include <sstream>

namespace pygcmc {
namespace system {
namespace molecular {

MolecularCombiner::MolecularCombiner() 
    : validator_(std::make_unique<MolecularValidator>()),
      matcher_(std::make_unique<MolecularMatcher>()),
      merger_(std::make_unique<MolecularMerger>()) {
}

std::shared_ptr<model::Molecular> MolecularCombiner::combine(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology) {
    
    // Validate input parameters
    validator_->validateCombination(structure, topology);

    // Create new Molecular object
    auto molecular = std::make_shared<model::Molecular>();

    // Copy data from Structure
    copyStructureData(molecular, structure);

    // Perform detailed matching
    matcher_->performDetailedMatching(molecular, topology);

    // Copy topology data
    copyTopologyData(molecular, topology);

    // Build lookup mappings
    buildLookupMappings(molecular, topology);

    return molecular;
}

std::shared_ptr<model::Molecular> MolecularCombiner::combineMultiple(
    const std::shared_ptr<model::Structure>& structure,
    const std::vector<std::shared_ptr<model::Topology>>& topologies) {
    
    if (!structure || topologies.empty()) {
        throw std::invalid_argument("Structure and Topologies cannot be null/empty");
    }

    // Create new Molecular object
    auto molecular = std::make_shared<model::Molecular>();

    // Copy data from Structure
    copyStructureData(molecular, structure);

    // Match residues with topologies - this is the key step that determines which topologies are actually used
    auto matched_topologies = matchResidues(structure, topologies);
    
    // Calculate totals AFTER matching
    size_t total_atoms = 0;
    size_t total_residues = 0;
    for (const auto& topology : matched_topologies) {
        total_atoms += topology->get_num_atoms();
        total_residues += topology->get_num_residues();
    }

    // Validate MATCHED topologies against structure
    validator_->validateMultipleCombination(structure, matched_topologies, total_atoms, total_residues);
    
    // Merge all matched topologies
    merger_->mergeTopologies(molecular, matched_topologies);
    
    return molecular;
}

void MolecularCombiner::copyStructureData(
    std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Structure>& structure) {
    
    molecular->atoms = structure->get_atoms();
    molecular->residues = structure->get_residues();
    molecular->terminals = structure->get_terminals();
    molecular->helices = structure->get_helices();
    molecular->sheets = structure->get_sheets();
    molecular->ssbonds = structure->get_ssbonds();
    molecular->boxDimensions = structure->get_box_dimensions();
}

void MolecularCombiner::copyTopologyData(
    std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Topology>& topology) {
    
    const auto num_atoms = static_cast<size_t>(topology->get_num_atoms());
    const auto num_residues = static_cast<size_t>(topology->get_num_residues());
    const auto num_segments = static_cast<size_t>(topology->get_num_segments());

    molecular->topology_atoms.reserve(num_atoms);
    molecular->topology_residues.reserve(num_residues);
    molecular->segments.reserve(num_segments);

    for (size_t i = 0; i < num_atoms; ++i) {
        molecular->topology_atoms.push_back(topology->get_atom(static_cast<int>(i)));
    }
    for (size_t i = 0; i < num_residues; ++i) {
        molecular->topology_residues.push_back(topology->get_residue(static_cast<int>(i)));
    }
    for (size_t i = 0; i < num_segments; ++i) {
        molecular->segments.push_back(topology->get_segment(static_cast<int>(i)));
    }

    // Copy bonding information
    molecular->bonds = topology->get_bonds();
    molecular->angles = topology->get_angles();
    molecular->dihedrals = topology->get_dihedrals();
    molecular->donors = topology->get_donors();
    molecular->acceptors = topology->get_acceptors();
    molecular->exclusions = topology->get_exclusions();
    molecular->groups = topology->get_groups();
    molecular->cmaps = topology->get_cmaps();
    
    // Add standardized CMAPs
    for (const auto& cmap : topology->get_cmaps()) {
        molecular->add_standard_cmap(cmap);
    }
    
    molecular->titles = topology->get_titles();
}

void MolecularCombiner::buildLookupMappings(
    std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Topology>& topology) {
    
    const auto num_atoms = static_cast<size_t>(topology->get_num_atoms());
    const auto num_residues = static_cast<size_t>(topology->get_num_residues());
    const auto num_segments = static_cast<size_t>(topology->get_num_segments());

    // Copy lookup mappings
    for (size_t i = 0; i < num_segments; ++i) {
        const auto& segment = topology->get_segment(static_cast<int>(i));
        molecular->segment_map[segment.name] = segment.id;
    }

    for (size_t i = 0; i < num_residues; ++i) {
        const auto& residue = topology->get_residue(static_cast<int>(i));
        molecular->residue_map[std::make_pair(residue.name, residue.number)] = residue.id;
    }

    for (size_t i = 0; i < num_atoms; ++i) {
        const auto& atom = topology->get_atom(static_cast<int>(i));
        const auto& residue = topology->get_residue(atom.residue_id);
        molecular->atom_map[std::make_tuple(residue.name, residue.number, atom.name)] = atom.id;
    }
}

std::vector<std::shared_ptr<model::Topology>> MolecularCombiner::matchResidues(
    const std::shared_ptr<model::Structure>& structure,
    const std::vector<std::shared_ptr<model::Topology>>& topologies) {
    
    const auto& residues = structure->get_residues();
    std::vector<bool> residue_matched(residues.size(), false);
    std::vector<std::shared_ptr<model::Topology>> matched_topologies;
    
    // Try to match each topology
    for (size_t start_idx = 0; start_idx < residues.size(); ++start_idx) {
        if (residue_matched[start_idx]) continue;
        
        for (const auto& topology : topologies) {
            size_t matched_count = 0;
            if (matcher_->matchResidueSequence(residues, topology, start_idx, matched_count)) {
                // Mark matched residues
                for (size_t i = 0; i < matched_count; ++i) {
                    residue_matched[start_idx + i] = true;
                }
                matched_topologies.push_back(topology);
                break;
            }
        }
    }
    
    // Check if all residues have been matched
    for (size_t i = 0; i < residue_matched.size(); ++i) {
        if (!residue_matched[i]) {
            std::stringstream ss;
            ss << "Could not find matching topology for residue " 
               << residues[i]->get_resname()
               << " " << residues[i]->get_ires();
            throw std::runtime_error(ss.str());
        }
    }
    
    return matched_topologies;
}

} // namespace molecular
} // namespace system
} // namespace pygcmc 