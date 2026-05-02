#include "SystemMain.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace system {

void System::initialize_parameters() {
    auto& file_info = params_.get_file_info();
    auto& fragment_info = params_.get_fragment_info();
    auto& basic_info = params_.get_basic_info();
    auto& energy_info = params_.get_energy_info();
    auto& mc_info = params_.get_mc_info();
    auto& bias_info = params_.get_bias_info();

    // Set water molecule index and density
    for (size_t i = 0; i < file_info.fragment_names.size(); i++) {
        if (file_info.fragment_names[i] == "sol") {
            fragment_info.water_density = fragment_info.conc_list[i];
            fragment_info.water_index = i;
            break;
        }
    }

    // Check box definition
    if (!basic_info.is_box) {
        throw std::runtime_error("Box size not defined!");
    }

    // Process cavity list
    process_cavity_list();

    // Initialize MC time list
    initialize_mc_time_list();

    // Calculate pair list cutoff
    energy_info.pairlist_cutoff = energy_info.fragment_cutoff +
        mc_info.max_translation_dist * std::sqrt(3.0f) + 3.0f;
    energy_info.pairlist_cutoff_squared =
        energy_info.pairlist_cutoff * energy_info.pairlist_cutoff;

    // Calculate beta value
    mc_info.beta = 1.0f / (mc_info.BOLTZMANN * mc_info.temperature);

    // Check configuration bias parameters
    if (bias_info.use_conf_bias && bias_info.num_conf_bias_trials < 1) {
        throw std::runtime_error("num_conf_bias_trial needs to be greater than 0");
    }
}

void System::process_cavity_list() {
    auto& fragment_info = params_.get_fragment_info();

    for (float radius : fragment_info.radius_list) {
        auto it = std::lower_bound(fragment_info.cavity_list.begin(),
                                 fragment_info.cavity_list.end(), radius);
        if (it == fragment_info.cavity_list.end() || *it != radius) {
            fragment_info.cavity_list.insert(it, radius);
        }
    }

    fragment_info.cavity_index_list.resize(fragment_info.radius_list.size());
    for (size_t i = 0; i < fragment_info.radius_list.size(); i++) {
        auto it = std::find(fragment_info.cavity_list.begin(),
                          fragment_info.cavity_list.end(),
                          fragment_info.radius_list[i]);
        fragment_info.cavity_index_list[i] =
            it != fragment_info.cavity_list.end() ?
            std::distance(fragment_info.cavity_list.begin(), it) : -1;
    }
}

void System::initialize_mc_time_list() {
    auto& file_info = params_.get_file_info();
    auto& mc_info = params_.get_mc_info();

    if (mc_info.mc_time_list.empty()) {
        mc_info.mc_time_list.resize(file_info.fragment_names.size(),
            1.0f / file_info.fragment_names.size());
    }

    mc_info.mc_time_cumulative.resize(file_info.fragment_names.size());
    mc_info.mc_time_cumulative[0] = mc_info.mc_time_list[0];
    for (size_t i = 1; i < file_info.fragment_names.size(); i++) {
        mc_info.mc_time_cumulative[i] =
            mc_info.mc_time_cumulative[i-1] + mc_info.mc_time_list[i];
    }

    // Normalize
    float total = mc_info.mc_time_cumulative.back();
    if (total != 1.0f) {
        for (float& time : mc_info.mc_time_cumulative) {
            time /= total;
        }
    }
}

} // namespace system
} // namespace pygcmc
