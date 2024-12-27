// modules/core/include/pygcmc/core/io/ff_parser.hpp

#ifndef PYGCMC_CORE_IO_FF_PARSER_HPP
#define PYGCMC_CORE_IO_FF_PARSER_HPP

#include "parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief Force Field Parser
 * 
 * 解析力场参数文件（如 TOP 文件中的非键参数）。
 */
class FFParser {
public:
    /**
     * @brief 解析力场文件
     * 
     * @param filename 力场文件路径
     * @return std::pair<NBMap, NBFixMap> 包含非键参数和修正参数的映射
     */
    static std::pair<NBMap, NBFixMap> parse(const std::string& filename);

private:
    /**
     * @brief 解析 [ atomtypes ] 部分
     * 
     * @param line 当前行内容
     * @param nb_dict 非键参数映射
     * @return true 解析成功
     * @return false 解析失败
     */
    static bool parse_atomtypes_line(const std::string& line, NBMap& nb_dict);

    /**
     * @brief 解析 [ nonbond_params ] 部分
     * 
     * @param line 当前行内容
     * @param nbfix_dict 修正参数映射
     * @return true 解析成功
     * @return false 解析失败
     */
    static bool parse_nonbond_params_line(const std::string& line, NBFixMap& nbfix_dict);

    /**
     * @brief 解析 [ pairtypes ] 部分
     * 
     * @param line 当前行内容
     * @param nbfix_dict 修正参数映射
     * @return true 解析成功
     * @return false 解析失败
     */
    static bool parse_pairtypes_line(const std::string& line, NBFixMap& nbfix_dict);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_FF_PARSER_HPP
