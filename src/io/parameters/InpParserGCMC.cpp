#include "InpParserGCMC.hpp"
#include <fstream>
#include <sstream>

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
}

void InpParserGCMC::parse_line_ext(const std::string& key, const std::string& value, model::param::Param& param) {
    auto& mc_info = param.get_mc_info();
    auto& frag_info = param.get_fragment_info();
    auto& bias_info = param.get_bias_info();
    auto& space_info = param.get_space_info();
    auto& energy_info = param.get_energy_info();

    if (key == "mctime") {
        // Support accumulation of multiple mctime lines
        auto times = InpParserStructures::parse_float_vector(value);
        for (float t : times) {
            mc_info.mc_time_list.push_back(t);
        }
    } else if (key == "fragradius") {
        frag_info.radius_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "fragconf" || key == "fragconfs") {
        // keep both for compatibility
        frag_info.conf_list = InpParserStructures::parse_int_vector(value);
        frag_info.fragconf_list = frag_info.conf_list;
    } else if (key == "num_conf_bias_trial" || key == "confbias_trials") {
        bias_info.num_conf_bias_trials = static_cast<unsigned int>(std::stoi(value));
    } else if (key == "cavity_grid_dx") {
        // Map to grid_dx if given (fallback)
        space_info.grid_spacing = std::stof(value);
    } else if (key == "probe_radius") {
        // Map to sigma (approximate) if present
        const float r = std::stof(value);
        bias_info.sigma = r;
        bias_info.sigma_squared = r * r;
    } else if (key == "wdens") {
        // Water density output control
        mc_info.wdens = std::stof(value);
    } else if (key == "target_numwaters" || key == "target_num_waters") {
        // Target number of water molecules - store in both places
        frag_info.target_num_waters = std::stoi(value);
    } else if (key == "gcmc_region") {
        // GCMC insertion region (sphere/box specification)
        space_info.gcmc_region = value;
    } else if (key == "exclude_protein_volume") {
        // Exclude protein volume from cavity bias
        space_info.exclude_protein_volume = (value == "yes" || value == "true" || value == "1");
    } else if (key == "use_vdw_radius_for_grid") {
        // Use VDW radii for grid generation
        space_info.use_vdw_radius_for_grid = (value == "yes" || value == "true" || value == "1");
    } else if (key == "exclude_hydrogens_from_grid") {
        // Exclude hydrogens from grid occupancy
        space_info.exclude_hydrogens_from_grid = (value == "yes" || value == "true" || value == "1");
    } else if (key == "use_switching") {
        // Enable switching function
        mc_info.use_switching = (value == "yes" || value == "true" || value == "1");
    } else if (key == "switch_r_on" || key == "switch_ron") {
        // Switching function r_on
        mc_info.switch_r_on = std::stof(value);
    } else if (key == "switch_r_off" || key == "switch_roff") {
        // Switching function r_off
        mc_info.switch_r_off = std::stof(value);
    } else if (key == "pairlist_freq") {
        // Pairlist update frequency
        energy_info.pairlist_freq = static_cast<unsigned int>(std::stoi(value));
    }
}

void InpParserGCMC::enhance_param(model::param::Param& param) {
    // Ensure mc_time_cumulative if mc_time_list provided
    auto& mc_info = param.get_mc_info();
    if (!mc_info.mc_time_list.empty()) {
        mc_info.mc_time_cumulative.clear();
        float s = 0.0f;
        for (float w : mc_info.mc_time_list) {
            s += w;
            mc_info.mc_time_cumulative.push_back(s);
        }
        if (s > 0.0f) {
            for (auto& v : mc_info.mc_time_cumulative) v /= s;
        }
    }
}

} // namespace parameters
} // namespace io
} // namespace pygcmc