#include "InpParserGCMC.hpp"
#include <algorithm>
#include <string>
#include <vector>

namespace pygcmc {
namespace io {
namespace parameters {

void InpParserGCMC::parse_line_ext(const std::string& key, const std::string& value, model::param::Param& param) {
    bool handled = parse_line_ext_core(key, value, param);
    if (!handled) {
        handled = parse_line_ext_legacy(key, value, param);
    }

    if (!handled) {
        return;
    }

    auto& basic_info = param.get_basic_info();
    auto& handled_keys = basic_info.inp_keys_handled;
    if (std::find(handled_keys.begin(), handled_keys.end(), key) == handled_keys.end()) {
        handled_keys.push_back(key);
    }
}

} // namespace parameters
} // namespace io
} // namespace pygcmc

