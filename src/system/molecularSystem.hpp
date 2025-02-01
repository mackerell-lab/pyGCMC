#pragma once
#ifndef PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP
#define PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP

#include <memory>
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

private:
    std::shared_ptr<model::Molecular> molecular_;
};

} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_MOLECULAR_SYSTEM_HPP
