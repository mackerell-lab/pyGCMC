#pragma once
#ifndef PYGCMC_MODEL_MOLECULAR_HPP
#define PYGCMC_MODEL_MOLECULAR_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_map>
#include <array>
#include "model/atom.hpp"
#include "model/residue.hpp"
#include "model/structure.hpp"
#include "model/topology.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief 标准化的CMAP结构，用于统一处理PSF和TOP格式
 */
struct StandardCmap {
    std::array<int, 5> atoms;     ///< 标准化的5个原子索引
    std::array<int, 8> raw_atoms; ///< 原始格式的原子索引（8个for PSF，5个+3个-1 for TOP）
    bool is_psf_format;           ///< 是否是PSF格式
    int function_type = 1;        ///< CMAP函数类型
};

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
    std::vector<TopologyCmap> cmaps;  // 原始CMAP项
    std::vector<StandardCmap> standard_cmaps;  // 标准化的CMAP项
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

    // 标准化CMAP相关方法
    void add_standard_cmap(const TopologyCmap& cmap) {
        StandardCmap std_cmap;
        std_cmap.raw_atoms = cmap.atoms;
        std_cmap.is_psf_format = (cmap.atoms[5] != -1);  // 判断是否为PSF格式
        
        // 设置标准化的5个原子
        if (std_cmap.is_psf_format) {
            // PSF格式：使用第1-4个原子和第8个原子
            for (int i = 0; i < 4; ++i) {
                std_cmap.atoms[i] = cmap.atoms[i];
            }
            std_cmap.atoms[4] = cmap.atoms[7];  // 使用第8个原子作为第5个原子
        } else {
            // TOP格式：直接使用前5个原子
            for (int i = 0; i < 5; ++i) {
                std_cmap.atoms[i] = cmap.atoms[i];
            }
        }
        std_cmap.function_type = cmap.function_type;
        standard_cmaps.push_back(std_cmap);
    }

    // 获取标准化CMAP数量
    size_t get_num_standard_cmaps() const { return standard_cmaps.size(); }

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
        standard_cmaps.clear();
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
