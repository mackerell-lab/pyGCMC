#include "system/molecularSystem.hpp"
#include <stdexcept>

namespace pygcmc {
namespace system {

std::shared_ptr<model::Molecular> MolecularSystem::combine(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology) {
    
    if (!structure || !topology) {
        throw std::invalid_argument("Structure and Topology pointers cannot be null");
    }

    // 创建新的Molecular对象
    molecular_ = std::make_shared<model::Molecular>();

    // 从Structure复制数据
    molecular_->atoms = structure->get_atoms();
    molecular_->residues = structure->get_residues();
    molecular_->terminals = structure->get_terminals();
    molecular_->helices = structure->get_helices();
    molecular_->sheets = structure->get_sheets();
    molecular_->ssbonds = structure->get_ssbonds();
    molecular_->boxDimensions = structure->get_box_dimensions();

    // 从Topology复制数据
    // 复制原子、残基和片段信息
    const auto num_atoms = static_cast<size_t>(topology->get_num_atoms());
    const auto num_residues = static_cast<size_t>(topology->get_num_residues());
    const auto num_segments = static_cast<size_t>(topology->get_num_segments());

    molecular_->topology_atoms.reserve(num_atoms);
    molecular_->topology_residues.reserve(num_residues);
    molecular_->segments.reserve(num_segments);

    for (size_t i = 0; i < num_atoms; ++i) {
        molecular_->topology_atoms.push_back(topology->get_atom(static_cast<int>(i)));
    }
    for (size_t i = 0; i < num_residues; ++i) {
        molecular_->topology_residues.push_back(topology->get_residue(static_cast<int>(i)));
    }
    for (size_t i = 0; i < num_segments; ++i) {
        molecular_->segments.push_back(topology->get_segment(static_cast<int>(i)));
    }

    // 复制键合信息
    molecular_->bonds = topology->get_bonds();
    molecular_->angles = topology->get_angles();
    molecular_->dihedrals = topology->get_dihedrals();
    molecular_->donors = topology->get_donors();
    molecular_->acceptors = topology->get_acceptors();
    molecular_->exclusions = topology->get_exclusions();
    molecular_->groups = topology->get_groups();
    molecular_->cmaps = topology->get_cmaps();
    molecular_->titles = topology->get_titles();

    // 复制查找映射
    for (size_t i = 0; i < num_segments; ++i) {
        const auto& segment = topology->get_segment(static_cast<int>(i));
        molecular_->segment_map[segment.name] = segment.id;
    }

    for (size_t i = 0; i < num_residues; ++i) {
        const auto& residue = topology->get_residue(static_cast<int>(i));
        molecular_->residue_map[std::make_pair(residue.name, residue.number)] = residue.id;
    }

    for (size_t i = 0; i < num_atoms; ++i) {
        const auto& atom = topology->get_atom(static_cast<int>(i));
        const auto& residue = topology->get_residue(atom.residue_id);
        molecular_->atom_map[std::make_tuple(residue.name, residue.number, atom.name)] = atom.id;
    }

    // 验证数据一致性
    if (molecular_->atoms.size() != num_atoms) {
        throw std::runtime_error("Inconsistent number of atoms between Structure and Topology");
    }

    if (molecular_->residues.size() != num_residues) {
        throw std::runtime_error("Inconsistent number of residues between Structure and Topology");
    }

    return molecular_;
}

} // namespace system
} // namespace pygcmc
