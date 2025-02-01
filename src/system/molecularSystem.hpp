#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP
#define PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP

#include <memory>
#include <vector>
#include "model/molecular.hpp"
#include "model/structure.hpp"
#include "model/topology.hpp"

namespace pygcmc {
namespace system {

/**
 * @brief 分子系统类，用于管理和操作分子数据
 */
class MolecularSystem {
public:
    MolecularSystem() = default;
    ~MolecularSystem() = default;

    /**
     * @brief 合并Structure和Topology数据到Molecular对象
     * @param structure 结构数据
     * @param topology 拓扑数据
     * @return 合并后的Molecular对象
     */
    std::shared_ptr<model::Molecular> combine(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology);

    /**
     * @brief 合并多个Structure和Topology数据到Molecular对象
     * @param structure 结构数据
     * @param topologies 拓扑数据
     * @return 合并后的Molecular对象
     */
    std::shared_ptr<model::Molecular> combine_multiple(
        const std::shared_ptr<model::Structure>& structure,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);

private:
    std::shared_ptr<model::Molecular> molecular_;

    /**
     * @brief 检查残基序列是否匹配
     * @param pdb_residues 结构中的残基序列
     * @param topology 拓扑数据
     * @param start_idx 起始索引
     * @param matched_count 匹配的残基数量
     * @return 是否匹配
     */
    bool match_residue_sequence(
        const std::vector<std::shared_ptr<model::Residue>>& pdb_residues,
        const std::shared_ptr<model::Topology>& topology,
        size_t start_idx,
        size_t& matched_count);

    /**
     * @brief 合并多个拓扑文件
     * @param molecular 分子对象
     * @param topologies 拓扑数据
     */
    void merge_topologies(
        std::shared_ptr<model::Molecular>& molecular,
        const std::vector<std::shared_ptr<model::Topology>>& topologies);

    /**
     * @brief 验证原子类型匹配
     * @param pdb_res 结构中的残基
     * @param top_res 拓扑中的残基
     * @param topology 拓扑数据
     */
    void verify_atom_types(
        const std::shared_ptr<model::Residue>& pdb_res,
        const model::TopologyResidue& top_res,
        const std::shared_ptr<model::Topology>& topology);
};

} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP
