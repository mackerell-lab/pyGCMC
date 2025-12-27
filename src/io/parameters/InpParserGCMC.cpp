#include "InpParserGCMC.hpp"
#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace pygcmc {
namespace io {
namespace parameters {

using ::pygcmc::io::parameters::InpParserStructures;

void InpParserGCMC::parse_to_param(const std::string& filename, model::param::Param& param) {
    // First parse with the base parser
    InpParserMain::parse_to_param(filename, param);

    // Then parse GCMC-specific keys
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open input file: " + filename);
    }
    std::string line;
    while (std::getline(file, line)) {
        line = InpParserStructures::trim(line);
        if (line.empty() || line[0] == '#') continue;
        auto tokens = InpParserStructures::split(line, ':');
        if (tokens.size() != 2) continue;
        const std::string key = InpParserStructures::trim(tokens[0]);
        const std::string value = InpParserStructures::trim(tokens[1]);
        parse_line_ext(key, value, param);
    }
    enhance_param(param);

    // Derive unknown INP keys after both base + GCMC-specific parsing passes.
    auto& basic_info = param.get_basic_info();
    basic_info.inp_keys_unknown.clear();
    for (const auto& k : basic_info.inp_keys_seen) {
        if (std::find(basic_info.inp_keys_handled.begin(), basic_info.inp_keys_handled.end(), k) ==
            basic_info.inp_keys_handled.end()) {
            basic_info.inp_keys_unknown.push_back(k);
        }
    }
}

} // namespace parameters
} // namespace io
} // namespace pygcmc
