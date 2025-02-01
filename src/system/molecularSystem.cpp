#include "system/molecularSystem.hpp"
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <map>

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

    // 验证原子总数
    if (molecular_->atoms.size() != num_atoms) {
        std::stringstream ss;
        ss << "Inconsistent total number of atoms: Structure has " 
           << molecular_->atoms.size() << " atoms, but Topology has " 
           << num_atoms << " atoms";
        throw std::runtime_error(ss.str());
    }

    // 验证每个残基的原子数和类型
    for (const auto& res : molecular_->residues) {
        const auto& res_atoms = res->get_atoms();
        bool found_matching_res = false;
        
        // 在topology中查找匹配的残基
        for (size_t i = 0; i < num_residues; ++i) {
            const auto& top_res = topology->get_residue(static_cast<int>(i));
            
            // 检查是否是蛋白质残基
            const std::string& resname = res->get_resname();
            bool is_protein = (resname.length() == 3) && 
                            (resname != "SOL") && (resname != "WAT") && (resname != "HOH") &&  // 水分子
                            (resname != "ION") && (resname != "CLA") && (resname != "SOD") &&  // 离子
                            (resname != "TIP") && (resname != "SPC");  // 其他水模型
            
            bool residue_matches;
            if (is_protein) {
                residue_matches = (top_res.name == resname && top_res.number == res->get_ires());
            } else {
                residue_matches = (top_res.name == resname);
            }
            
            if (residue_matches) {
                found_matching_res = true;
                
                // 检查原子数量
                if (res_atoms.size() != top_res.atoms.size()) {
                    std::stringstream ss;
                    ss << "Inconsistent number of atoms in residue " << res->get_resname() 
                       << " " << res->get_ires() << ": Structure has " 
                       << res_atoms.size() << " atoms, but Topology has " 
                       << top_res.atoms.size() << " atoms";
                    throw std::runtime_error(ss.str());
                }

                // 检查原子类型
                for (size_t j = 0; j < res_atoms.size(); ++j) {
                    const auto& pdb_atom = res_atoms[j];
                    const auto& top_atom = topology->get_atom(top_res.atoms[j]);
                    
                    // 从PDB原子名称中提取元素
                    std::string pdb_element = pdb_atom->get_element();
                    if (pdb_element.empty()) {
                        pdb_element = pdb_atom->get_type();
                    }
                    
                    // 从topology原子类型中提取元素
                    std::string top_element = top_atom.type;
                    
                    // 比较第一个字母（转为大写）
                    char pdb_first = std::toupper(pdb_element[0]);
                    char top_first = std::toupper(top_element[0]);
                    
                    if (pdb_first != top_first) {
                        std::stringstream ss;
                        ss << "Mismatched atom elements in residue " << res->get_resname() 
                           << " " << res->get_ires() << ": Structure has " 
                           << pdb_first << " (from " << pdb_atom->get_type() 
                           << "), but Topology has " << top_first 
                           << " (from " << top_atom.type << ")";
                        throw std::runtime_error(ss.str());
                    }
                }
                break;
            }
        }
        
        if (!found_matching_res) {
            std::stringstream ss;
            ss << "Could not find matching residue in topology for " 
               << res->get_resname() << " " << res->get_ires();
            throw std::runtime_error(ss.str());
        }
    }

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

    return molecular_;
}

} // namespace system
} // namespace pygcmc
