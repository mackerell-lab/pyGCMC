#include "MolecularMerger.hpp"
#include <stdexcept>
#include <sstream>

namespace pygcmc {
namespace system {
namespace molecular {

void MolecularMerger::mergeTopologies(
    std::shared_ptr<model::Molecular>& molecular,
    const std::vector<std::shared_ptr<model::Topology>>& topologies) {
    
    // Clear existing topology data
    molecular->topology_atoms.clear();
    molecular->topology_residues.clear();
    molecular->segments.clear();
    molecular->bonds.clear();
    molecular->angles.clear();
    molecular->dihedrals.clear();
    molecular->donors.clear();
    molecular->acceptors.clear();
    molecular->exclusions.clear();
    molecular->groups.clear();
    
    // Save existing CMAP
    std::vector<model::TopologyCmap> existing_cmaps = molecular->cmaps;
    molecular->cmaps.clear();
    
    size_t atom_offset = 0;
    size_t residue_offset = 0;
    
    for (const auto& topology : topologies) {
        // Copy atoms and residues
        copyAtomsWithOffset(molecular, topology, atom_offset, residue_offset);
        copyResiduesWithOffset(molecular, topology, atom_offset, residue_offset);
        
        // Copy bonding information
        copyBondingWithOffset(molecular, topology, atom_offset);
        
        // Copy CMAP information
        copyCMAPWithOffset(molecular, topology, atom_offset, existing_cmaps);
        
        // Update offsets
        atom_offset += topology->get_num_atoms();
        residue_offset += topology->get_num_residues();
    }
}

void MolecularMerger::copyAtomsWithOffset(
    std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Topology>& topology,
    size_t atom_offset,
    size_t residue_offset) {
    
    for (int i = 0; i < topology->get_num_atoms(); ++i) {
        auto atom = topology->get_atom(i);
        atom.id += atom_offset;
        atom.residue_id += residue_offset;
        molecular->topology_atoms.push_back(atom);
    }
}

void MolecularMerger::copyResiduesWithOffset(
    std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Topology>& topology,
    size_t atom_offset,
    size_t residue_offset) {
    
    for (int i = 0; i < topology->get_num_residues(); ++i) {
        auto residue = topology->get_residue(i);
        residue.id += residue_offset;
        for (auto& atom_id : residue.atoms) {
            atom_id += atom_offset;
        }
        molecular->topology_residues.push_back(residue);
    }
}

void MolecularMerger::copyBondingWithOffset(
    std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Topology>& topology,
    size_t atom_offset) {
    
    // Copy bond information
    for (const auto& bond : topology->get_bonds()) {
        model::TopologyBond new_bond = bond;
        new_bond.atom1 += atom_offset;
        new_bond.atom2 += atom_offset;
        molecular->bonds.push_back(new_bond);
    }
    
    // Copy angle information
    for (const auto& angle : topology->get_angles()) {
        model::TopologyAngle new_angle = angle;
        new_angle.atom1 += atom_offset;
        new_angle.atom2 += atom_offset;
        new_angle.atom3 += atom_offset;
        molecular->angles.push_back(new_angle);
    }
    
    // Copy dihedral information
    for (const auto& dihedral : topology->get_dihedrals()) {
        model::TopologyDihedral new_dihedral = dihedral;
        new_dihedral.atom1 += atom_offset;
        new_dihedral.atom2 += atom_offset;
        new_dihedral.atom3 += atom_offset;
        new_dihedral.atom4 += atom_offset;
        molecular->dihedrals.push_back(new_dihedral);
    }
}

void MolecularMerger::copyCMAPWithOffset(
    std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Topology>& topology,
    size_t atom_offset,
    const std::vector<model::TopologyCmap>& existing_cmaps) {
    
    // Copy CMAP information
    for (const auto& cmap : topology->get_cmaps()) {
        model::TopologyCmap new_cmap = cmap;
        
        // Update indices for all 8 atoms
        for (size_t i = 0; i < new_cmap.atoms.size(); ++i) {
            if (new_cmap.atoms[i] >= 0) {
                const auto& atom = topology->get_atom(new_cmap.atoms[i]);
                const auto& res = topology->get_residue(atom.residue_id);
                if (res.name != "SOL") {
                    new_cmap.atoms[i] += atom_offset;
                }
            }
        }
        
        molecular->cmaps.push_back(new_cmap);
        // Add standardized CMAP
        molecular->add_standard_cmap(new_cmap);
    }
    
    // Re-add previously saved CMAPs
    for (const auto& cmap : existing_cmaps) {
        molecular->cmaps.push_back(cmap);
        molecular->add_standard_cmap(cmap);
    }
}

} // namespace molecular
} // namespace system
} // namespace pygcmc 