// modules/core/include/pygcmc/core/utils.hpp

#ifndef PYGCMC_CORE_UTILS_HPP
#define PYGCMC_CORE_UTILS_HPP

#include <string>
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace core {
namespace utils {

/**
 * @brief 去除字符串首尾的空白字符
 * 
 * @param s 输入字符串
 * @return std::string 去除空白后的字符串
 */
inline std::string trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) {
        start++;
    }
    auto end = s.end();
    if (start != end) {
        do {
            end--;
        } while (std::distance(start, end) > 0 && std::isspace(*end));
        end++;
    }
    return std::string(start, end);
}

} // namespace utils
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_UTILS_HPP
