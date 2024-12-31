// modules/core/include/pygcmc/core/io/ff_parser.hpp

#ifndef PYGCMC_CORE_IO_FF_PARSER_HPP
#define PYGCMC_CORE_IO_FF_PARSER_HPP

#include <string>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <vector>
#include "pygcmc/core/io/parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

class FFParser {
public:
    // Parse CHARMM parameter file
    bool parse(const std::string& filename);

    // Get nonbonded parameters
    const std::unordered_map<std::string, ForceFieldPair>& get_nonbonded_params() const {
        return nonbonded_params_;
    }

    // Get NBFIX parameters
    const std::unordered_map<
        std::pair<std::string, std::string>, 
        ForceFieldPair, 
        PairStringHash
    >& get_nbfix_params() const {
        return nbfix_params_;
    }

    /**
     * @brief Update PDB atoms with force field parameters
     * @param atoms Vector of PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom>& atoms) const;

    /**
     * @brief Update PDB atoms with force field parameters
     * @param atoms Vector of PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom*>& atoms) const;   

    // 获取全局非键参数的访问器
    double get_cutnb() const { return cutnb_; }
    double get_ctofnb() const { return ctofnb_; }
    double get_ctonnb() const { return ctonnb_; }
    double get_eps() const { return eps_; }
    double get_e14fac() const { return e14fac_; }
    double get_wmin() const { return wmin_; }

private:
    // Helper functions for parsing sections
    void parse_nonbonded_line(const std::string& line);
    void parse_nbfix_line(const std::string& line);

    // 存储原子类型的非键参数
    std::unordered_map<std::string, ForceFieldPair> nonbonded_params_;

    // 存储原子类型对的 NBFIX 参数
    std::unordered_map<
        std::pair<std::string, std::string>, 
        ForceFieldPair, 
        PairStringHash
    > nbfix_params_;

    // 全局非键参数，设置默认值
    double cutnb_ = 14.0;    // 非键相互作用截断距离
    double ctofnb_ = 12.0;   // 开始切换的距离
    double ctonnb_ = 10.0;   // 开始计算的距离
    double eps_ = 1.0;       // 介电常数
    double e14fac_ = 1.0;    // 1-4相互作用的缩放因子
    double wmin_ = 1.5;      // 最小权重
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_FF_PARSER_HPP

