#pragma once
#ifndef PYGCMC_IO_PARAMETERS_INPPARSERGCMC_HPP
#define PYGCMC_IO_PARAMETERS_INPPARSERGCMC_HPP

#include "InpParserMain.hpp"
#include "InpParserStructures.hpp"
#include "../../model/param/ParamMain.hpp"
#include <string>

namespace pygcmc {
namespace io {
namespace parameters {

class InpParserGCMC {
public:
    // Parse base keys with InpParserMain, then parse GCMC-specific keys
    static void parse_to_param(const std::string& filename, model::param::Param& param);

    // Optional post-processing (derive/validate after all keys present)
    static void enhance_param(model::param::Param& param);

private:
    static void parse_line_ext(const std::string& key, const std::string& value, model::param::Param& param);

    // Split large parsing/normalization logic into smaller translation units.
    static bool parse_line_ext_core(const std::string& key, const std::string& value, model::param::Param& param);
    static bool parse_line_ext_legacy(const std::string& key, const std::string& value, model::param::Param& param);

    // `enhance_param` implementation helpers.
    // Returns true when the unit mode is treated as native nm/kJ (used for heuristic warnings).
    static bool enhance_param_units(model::param::Param& param);
    static void enhance_param_warnings(model::param::Param& param, bool nm_mode);
    static void enhance_param_finalize(model::param::Param& param);
};

} // namespace parameters
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_PARAMETERS_INPPARSERGCMC_HPP
