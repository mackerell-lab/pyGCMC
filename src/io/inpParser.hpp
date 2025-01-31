#pragma once
#ifndef PYGCMC_IO_INPPARSER_HPP
#define PYGCMC_IO_INPPARSER_HPP

#include <string>
#include "model/param.hpp"

namespace pygcmc {
namespace io {

/**
 * @brief GCMC输入文件解析器类
 * @details 用于解析GCMC模拟的输入文件，将参数存储到Param数据结构中
 * 支持的参数包括：
 * - 文件路径参数：par, fragitp, atomtypes, top, pdb等
 * - 空间参数：grid_dx, box_size, cutoff等
 * - 片段参数：fragname, fragconc, fragmuex
 * - 模拟控制参数：nprint, mcsteps等
 * - 偏置采样参数：use_cavity_bias, use_conf_bias
 */
class InpParser {
public:
    /**
     * @brief 解析输入文件并返回新的Param对象
     * @param filename 输入文件路径
     * @return model::Param 包含解析结果的Param对象
     * @throw std::runtime_error 如果文件不存在或解析失败
     */
    static model::Param parse_file(const std::string& filename);

    /**
     * @brief 解析输入字符串并返回新的Param对象
     * @param content 输入文件内容字符串
     * @return model::Param 包含解析结果的Param对象
     * @throw std::runtime_error 如果解析失败
     */
    static model::Param parse_string(const std::string& content);

    /**
     * @brief 解析输入文件并将结果存储到现有Param对象中
     * @param filename 输入文件路径
     * @param param 用于存储结果的Param对象
     * @throw std::runtime_error 如果文件不存在或解析失败
     */
    static void parse_to_param(const std::string& filename, model::Param& param);

    /**
     * @brief 解析输入字符串并将结果存储到现有Param对象中
     * @param content 输入文件内容字符串
     * @param param 用于存储结果的Param对象
     * @throw std::runtime_error 如果解析失败
     */
    static void parse_string_to_param(const std::string& content, model::Param& param);

private:
    /**
     * @brief 解析单行参数
     * @param key 参数名
     * @param value 参数值
     * @param param 用于存储结果的Param对象
     * @throw std::runtime_error 如果解析失败
     */
    static void parse_line(const std::string& key, const std::string& value, model::Param& param);

    /**
     * @brief 验证参数的有效性和一致性
     * @param param 需要验证的Param对象
     * @throw std::runtime_error 如果验证失败
     */
    static void validate_parameters(model::Param& param);
};

} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_INPPARSER_HPP
