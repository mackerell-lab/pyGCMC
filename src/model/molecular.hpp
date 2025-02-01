#pragma once
#ifndef PYGCMC_MODEL_MOLECULAR_HPP
#define PYGCMC_MODEL_MOLECULAR_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_map>
#include "model/atom.hpp"
#include "model/residue.hpp"
#include "model/structure.hpp"
#include "model/topology.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief 分子类，合并Structure和Topology的数据用于后续计算
 */
class Molecular {
public:
    Molecular() = default;
    ~Molecular() = default;

    // Structure (PDB) 相关数据
    std::vector<std::shared_ptr<Atom>> atoms;  // 原子列表
    std::vector<std::shared_ptr<Residue>> residues;  // 残基列表
    std::vector<Structure::TerminalInfo> terminals;  // 链终止信息
    std::map<std::string, std::vector<Structure::SecondaryStructure>> helices;  // 螺旋结构
    std::map<std::string, std::vector<std::string>> sheets;  // β片层结构
    std::vector<std::string> ssbonds;  // 二硫键
    std::vector<double> boxDimensions;  // 盒子尺寸

    // Topology (PSF/TOP) 相关数据
    std::vector<TopologyAtom> topology_atoms;  // 拓扑中的原子信息（包含电荷、质量等）
    std::vector<TopologyResidue> topology_residues;  // 拓扑中的残基信息
    std::vector<TopologySegment> segments;  // 片段信息
    std::vector<TopologyBond> bonds;  // 键
    std::vector<TopologyAngle> angles;  // 键角
    std::vector<TopologyDihedral> dihedrals;  // 二面角（包括improper）
    std::vector<TopologyDonor> donors;  // 氢键供体
    std::vector<TopologyAcceptor> acceptors;  // 氢键受体
    std::map<int, std::set<int>> exclusions;  // 非键排除
    std::vector<TopologyGroup> groups;  // 原子组
    std::vector<TopologyCmap> cmaps;  // CMAP项
    std::vector<std::string> titles;  // PSF文件的标题信息

    // 查找映射
    std::unordered_map<std::string, int> segment_map;  // segment_name -> index
    std::map<std::pair<std::string, int>, int> residue_map;  // (residue_name, number) -> index
    std::map<std::tuple<std::string, int, std::string>, int> atom_map;  // (residue_name, number, atom_name) -> index

    // 获取原子数量
    size_t get_num_atoms() const { return atoms.size(); }
    
    // 获取残基数量
    size_t get_num_residues() const { return residues.size(); }
    
    // 获取片段数量
    size_t get_num_segments() const { return segments.size(); }
    
    // 获取键数量
    size_t get_num_bonds() const { return bonds.size(); }
    
    // 获取键角数量
    size_t get_num_angles() const { return angles.size(); }
    
    // 获取二面角数量（不包括improper）
    size_t get_num_dihedrals() const {
        size_t count = 0;
        for (const auto& dihedral : dihedrals) {
            if (!dihedral.improper) count++;
        }
        return count;
    }
    
    // 获取improper数量
    size_t get_num_impropers() const {
        size_t count = 0;
        for (const auto& dihedral : dihedrals) {
            if (dihedral.improper) count++;
        }
        return count;
    }

    // 清除所有数据
    void clear() {
        // Structure数据
        atoms.clear();
        residues.clear();
        terminals.clear();
        helices.clear();
        sheets.clear();
        ssbonds.clear();
        boxDimensions.clear();

        // Topology数据
        topology_atoms.clear();
        topology_residues.clear();
        segments.clear();
        bonds.clear();
        angles.clear();
        dihedrals.clear();
        donors.clear();
        acceptors.clear();
        exclusions.clear();
        groups.clear();
        cmaps.clear();
        titles.clear();

        // 查找映射
        segment_map.clear();
        residue_map.clear();
        atom_map.clear();
    }
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULAR_HPP
