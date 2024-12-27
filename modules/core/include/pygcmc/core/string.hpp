// modules/core/include/pygcmc/core/string.hpp

#ifndef PYGCMC_CORE_UTILS_STRING_HPP
#define PYGCMC_CORE_UTILS_STRING_HPP

#include <string>
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace core {
namespace utils {
namespace string {

/**
 * @brief 将字符串转换为小写
 * 
 * @param s 输入字符串
 * @return std::string 小写后的字符串
 */
inline std::string to_lower(const std::string& s) {
    std::string result = s;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

/**
 * @brief 将字符串转换为大写
 * 
 * @param s 输入字符串
 * @return std::string 大写后的字符串
 */
inline std::string to_upper(const std::string& s) {
    std::string result = s;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return result;
}

// 添加更多字符串处理函数如有需要

} // namespace string
} // namespace utils
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_UTILS_STRING_HPP
