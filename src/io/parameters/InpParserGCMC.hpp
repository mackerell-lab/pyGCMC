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
};

} // namespace parameters
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_PARAMETERS_INPPARSERGCMC_HPP