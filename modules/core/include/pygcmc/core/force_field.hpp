// modules/core/include/pygcmc/core/force_field.hpp

#ifndef PYGCMC_CORE_FORCE_FIELD_HPP
#define PYGCMC_CORE_FORCE_FIELD_HPP

#include <string>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

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

namespace io {
    struct ForceFieldPair {
        double sigma;   ///< Lennard-Jones potential sigma parameter
        double epsilon; ///< Lennard-Jones potential epsilon parameter

        // Constructor with validation
        ForceFieldPair(double s = 0.0, double e = 0.0) 
            : sigma(s), epsilon(e) {
            if (!is_valid()) {
                throw ForceFieldError("Invalid force field parameters");
            }
        }

        // Validation method
        bool is_valid() const {
            return sigma > 0.0 && epsilon >= 0.0 &&
                   std::isfinite(sigma) && std::isfinite(epsilon);
        }

        /**
         * @brief Compute Lennard-Jones potential energy
         * @param distance Interatomic distance
         * @return Potential energy
         */
        double compute_lj_energy(double distance) const {
            if (distance < 1e-10) {
                return std::numeric_limits<double>::infinity();
            }
            double inv_r = sigma / distance;
            double inv_r6 = std::pow(inv_r, 6);
            double inv_r12 = inv_r6 * inv_r6;
            return 4.0 * epsilon * (inv_r12 - inv_r6);
        }

        /**
         * @brief Apply Lorentz-Berthelot combining rules
         */
        static ForceFieldPair lorentz_berthelot(const ForceFieldPair& a, const ForceFieldPair& b) {
            return ForceFieldPair(
                (a.sigma + b.sigma) * 0.5,
                std::sqrt(a.epsilon * b.epsilon)
            );
        }

        /**
         * @brief Apply geometric mean combining rules
         */
        static ForceFieldPair geometric_mean(const ForceFieldPair& a, const ForceFieldPair& b) {
            return ForceFieldPair(
                std::sqrt(a.sigma * b.sigma),
                std::sqrt(a.epsilon * b.epsilon)
            );
        }
    };
} // namespace io

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