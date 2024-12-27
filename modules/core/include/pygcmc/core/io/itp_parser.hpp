// modules/core/include/pygcmc/core/io/itp_parser.hpp

#ifndef PYGCMC_CORE_IO_ITP_PARSER_HPP
#define PYGCMC_CORE_IO_ITP_PARSER_HPP

#include "parser_common.hpp"
#include <vector>

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief ITP Parser
 * 
 * 解析 ITP 文件，获取分子信息。
 */
class ITPParser {
public:
    /**
     * @brief 解析 ITP 文件
     * 
     * @param filename ITP 文件路径
     * @return std::vector<ITPAtom> 解析得到的原子信息
     */
    static std::vector<ITPAtom> parse(const std::string& filename);

private:
    /**
     * @brief 解析 [ atoms ] 部分
     * 
     * @param is 输入流
     * @param atoms 存储解析得到的原子信息
     * @return true 解析成功
     * @return false 解析失败
     */
    static bool parse_atoms_section(std::istream& is, std::vector<ITPAtom>& atoms);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_ITP_PARSER_HPP
