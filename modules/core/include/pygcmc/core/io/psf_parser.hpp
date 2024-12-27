// modules/core/include/pygcmc/core/io/psf_parser.hpp

#ifndef PYGCMC_CORE_IO_PSF_PARSER_HPP
#define PYGCMC_CORE_IO_PSF_PARSER_HPP

#include "parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief PSF Parser
 * 
 * 解析 PSF 文件，获取系统的拓扑信息（如键）。
 */
class PSFParser {
public:
    /**
     * @brief 解析 PSF 文件
     * 
     * @param filename PSF 文件路径
     * @return PSFTopology 解析得到的拓扑信息
     */
    static PSFTopology parse(const std::string& filename);

private:
    /**
     * @brief 解析 BONDS 部分
     * 
     * @param is 输入流
     * @param topology 存储解析得到的拓扑信息
     * @return true 解析成功
     * @return false 解析失败
     */
    static bool parse_bonds_section(std::istream& is, PSFTopology& topology);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PSF_PARSER_HPP
