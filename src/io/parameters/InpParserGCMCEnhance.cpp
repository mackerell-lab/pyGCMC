#include "InpParserGCMC.hpp"

namespace pygcmc {
namespace io {
namespace parameters {

void InpParserGCMC::enhance_param(model::param::Param& param) {
    auto& basic_info = param.get_basic_info();
    if (!basic_info.inp_units_converted) {
        const bool nm_mode = enhance_param_units(param);
        enhance_param_warnings(param, nm_mode);
        basic_info.inp_units_converted = true;
    }
    enhance_param_finalize(param);
}

} // namespace parameters
} // namespace io
} // namespace pygcmc

