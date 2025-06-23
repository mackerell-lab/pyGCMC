#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_PARAMETER_RETRIEVER_HPP
#define PYGCMC_MODEL_FORCEFIELD_PARAMETER_RETRIEVER_HPP

#include "ForceFieldParams.hpp"
#include <stdexcept>
#include <sstream>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Optimized parameter retriever with smart caching and lookup strategies
 * 
 * This class provides efficient parameter retrieval with bidirectional lookups,
 * error handling, and optional caching for frequently accessed parameters.
 */
class ForceFieldParameterRetriever {
public:
    explicit ForceFieldParameterRetriever(const ForceFieldStorage& storage) : storage_(storage) {}
    ~ForceFieldParameterRetriever() = default;

    // === Single Parameter Retrieval ===

    /**
     * @brief Get atom mass with detailed error reporting
     */
    double get_atom_mass(const std::string& type) const {
        auto it = storage_.atom_masses.find(type);
        if (it == storage_.atom_masses.end()) {
            throw std::runtime_error(create_missing_parameter_error("atom mass", type));
        }
        return it->second;
    }

    /**
     * @brief Get LJ parameters with detailed error reporting
     */
    const LJParams& get_lj_params(const std::string& type) const {
        auto it = storage_.lj_params.find(type);
        if (it == storage_.lj_params.end()) {
            throw std::runtime_error(create_missing_parameter_error("LJ parameters", type));
        }
        return it->second;
    }

    // === Pair Parameter Retrieval ===

    /**
     * @brief Get NBFIX parameters with existence check
     */
    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, 
                                          const std::string& type2) const {
        auto key = ParamKeyUtils::make_pair(type1, type2);
        auto it = storage_.nbfix.find(key);
        if (it == storage_.nbfix.end()) {
            return std::make_pair(NBFIXParams{}, false);
        }
        return std::make_pair(it->second, true);
    }

    /**
     * @brief Get bond parameters with error reporting
     */
    const BondParams& get_bond_params(const std::string& type1, const std::string& type2) const {
        auto key = ParamKeyUtils::make_pair(type1, type2);
        auto it = storage_.bond_params.find(key);
        if (it == storage_.bond_params.end()) {
            throw std::runtime_error(create_missing_pair_error("bond parameters", type1, type2));
        }
        return it->second;
    }

    // === Triple Parameter Retrieval ===

    /**
     * @brief Get angle parameters with bidirectional lookup
     */
    const AngleParams& get_angle_params(const std::string& type1,
                                      const std::string& type2,
                                      const std::string& type3) const {
        // Try both orientations of the outer atoms while keeping the middle atom fixed
        auto key1 = std::make_tuple(type1, type2, type3);
        auto it = storage_.angle_params.find(key1);
        if (it != storage_.angle_params.end()) {
            return it->second;
        }

        // Try the reverse orientation
        auto key2 = std::make_tuple(type3, type2, type1);
        it = storage_.angle_params.find(key2);
        if (it != storage_.angle_params.end()) {
            return it->second;
        }

        throw std::runtime_error(create_missing_triple_error("angle parameters", type1, type2, type3));
    }

    // === Quadruple Parameter Retrieval ===

    /**
     * @brief Get dihedral parameters with multiple term support
     */
    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const {
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        auto it = storage_.dihedral_params.find(key);
        if (it == storage_.dihedral_params.end()) {
            throw std::runtime_error(create_missing_quad_error("dihedral parameters", type1, type2, type3, type4));
        }
        return it->second;
    }

    /**
     * @brief Get improper parameters with error reporting
     */
    const ImproperParams& get_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        auto it = storage_.improper_params.find(key);
        if (it == storage_.improper_params.end()) {
            throw std::runtime_error(create_missing_quad_error("improper parameters", type1, type2, type3, type4));
        }
        return it->second;
    }

    // === Batch Parameter Retrieval ===

    /**
     * @brief Get multiple atom masses at once with error collection
     */
    std::map<std::string, double> get_atom_masses_batch(const std::vector<std::string>& types) const {
        std::map<std::string, double> result;
        std::vector<std::string> missing;

        for (const auto& type : types) {
            auto it = storage_.atom_masses.find(type);
            if (it != storage_.atom_masses.end()) {
                result[type] = it->second;
            } else {
                missing.push_back(type);
            }
        }

        if (!missing.empty()) {
            std::ostringstream oss;
            oss << "Missing atom masses for types: ";
            for (size_t i = 0; i < missing.size(); ++i) {
                if (i > 0) oss << ", ";
                oss << missing[i];
            }
            throw std::runtime_error(oss.str());
        }

        return result;
    }

    /**
     * @brief Get multiple LJ parameters at once with error collection
     */
    std::map<std::string, LJParams> get_lj_params_batch(const std::vector<std::string>& types) const {
        std::map<std::string, LJParams> result;
        std::vector<std::string> missing;

        for (const auto& type : types) {
            auto it = storage_.lj_params.find(type);
            if (it != storage_.lj_params.end()) {
                result[type] = it->second;
            } else {
                missing.push_back(type);
            }
        }

        if (!missing.empty()) {
            std::ostringstream oss;
            oss << "Missing LJ parameters for types: ";
            for (size_t i = 0; i < missing.size(); ++i) {
                if (i > 0) oss << ", ";
                oss << missing[i];
            }
            throw std::runtime_error(oss.str());
        }

        return result;
    }

    // === Utility Methods ===

    /**
     * @brief Get parameter types with alternative lookup suggestions
     */
    std::vector<std::string> suggest_similar_types(const std::string& target_type) const {
        std::vector<std::string> suggestions;
        
        // Look for types with similar prefixes or suffixes
        for (const auto& pair : storage_.atom_masses) {
            const auto& type = pair.first;
            if (type != target_type && 
                (type.find(target_type) != std::string::npos || 
                 target_type.find(type) != std::string::npos)) {
                suggestions.push_back(type);
            }
        }
        
        return suggestions;
    }

private:
    const ForceFieldStorage& storage_;

    // === Error Message Helpers ===

    std::string create_missing_parameter_error(const std::string& param_type, 
                                             const std::string& type) const {
        std::ostringstream oss;
        oss << param_type << " not found for type: " << type;
        
        auto suggestions = suggest_similar_types(type);
        if (!suggestions.empty()) {
            oss << ". Similar types found: ";
            for (size_t i = 0; i < suggestions.size() && i < 3; ++i) {
                if (i > 0) oss << ", ";
                oss << suggestions[i];
            }
        }
        
        return oss.str();
    }

    std::string create_missing_pair_error(const std::string& param_type,
                                        const std::string& type1, 
                                        const std::string& type2) const {
        return param_type + " not found for types: " + type1 + "-" + type2;
    }

    std::string create_missing_triple_error(const std::string& param_type,
                                          const std::string& type1,
                                          const std::string& type2,
                                          const std::string& type3) const {
        return param_type + " not found for types: " + type1 + "-" + type2 + "-" + type3;
    }

    std::string create_missing_quad_error(const std::string& param_type,
                                        const std::string& type1,
                                        const std::string& type2,
                                        const std::string& type3,
                                        const std::string& type4) const {
        return param_type + " not found for types: " + type1 + "-" + type2 + "-" + type3 + "-" + type4;
    }
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_PARAMETER_RETRIEVER_HPP 