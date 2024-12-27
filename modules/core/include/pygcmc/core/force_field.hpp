// modules/core/include/pygcmc/core/force_field.hpp

#ifndef PYGCMC_CORE_FORCE_FIELD_HPP
#define PYGCMC_CORE_FORCE_FIELD_HPP

#include "pygcmc/core/io/parser_common.hpp"
#include <string>
#include <utility>
#include <unordered_map>
#include <stdexcept>

namespace pygcmc {
namespace core {

/**
 * @brief Force field related exceptions
 */
class ForceFieldError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/**
 * @brief Hash function for std::pair<std::string, std::string> keys in unordered_map
 */
struct PairStringHash {
    std::size_t operator()(const std::pair<std::string, std::string>& p) const {
        auto h1 = std::hash<std::string>{}(p.first);
        auto h2 = std::hash<std::string>{}(p.second);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

using NBMap = std::unordered_map<std::string, io::ForceFieldPair>;
using NBFixMap = std::unordered_map<std::pair<std::string, std::string>, 
                                   io::ForceFieldPair,
                                   PairStringHash>;

/**
 * @brief Validate non-bonded parameters
 */
inline bool validate_force_field(const NBMap& nb_params) {
    return std::all_of(nb_params.begin(), nb_params.end(),
                       [](const auto& pair) {
                           return pair.second.is_valid();
                       });
}

/**
 * @brief Validate non-bonded fix parameters
 */
inline bool validate_force_field(const NBFixMap& nbfix_params) {
    return std::all_of(nbfix_params.begin(), nbfix_params.end(),
                       [](const auto& pair) {
                           return pair.second.is_valid();
                       });
}

} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_FORCE_FIELD_HPP