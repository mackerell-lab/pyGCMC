#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_ANALYSIS_OPS_HPP
#define PYGCMC_MODEL_FORCEFIELD_ANALYSIS_OPS_HPP

#include "ForceFieldCore.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

// === Validation Method Implementations ===

inline bool ForceField::validate_force_field() const {
    return validator_.validate_force_field();
}

inline CompletenessResult ForceField::check_completeness(const std::set<std::string>& atom_types) const {
    return validator_.check_completeness(atom_types);
}

inline std::vector<std::string> ForceField::validate_consistency() const {
    return validator_.validate_consistency();
}

// === Analysis Method Implementations ===

inline ForceFieldStats ForceField::get_statistics() const {
    return analyzer_.get_statistics();
}

inline std::string ForceField::get_summary() const {
    return analyzer_.get_summary();
}

inline std::string ForceField::get_detailed_statistics() const {
    return analyzer_.get_detailed_statistics();
}

inline std::set<std::string> ForceField::find_missing_lj_params() const {
    return analyzer_.find_missing_lj_params();
}

inline std::set<std::string> ForceField::find_missing_masses() const {
    return analyzer_.find_missing_masses();
}

inline std::set<std::string> ForceField::get_all_referenced_types() const {
    return analyzer_.get_all_referenced_types();
}

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_ANALYSIS_OPS_HPP 