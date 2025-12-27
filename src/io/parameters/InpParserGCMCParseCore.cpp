#include "InpParserGCMC.hpp"
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace pygcmc {
namespace io {
namespace parameters {

using ::pygcmc::io::parameters::InpParserStructures;

bool InpParserGCMC::parse_line_ext_core(const std::string& key, const std::string& value, model::param::Param& param) {
    auto& mc_info = param.get_mc_info();
    auto& frag_info = param.get_fragment_info();
    auto& bias_info = param.get_bias_info();
    auto& space_info = param.get_space_info();
    auto& energy_info = param.get_energy_info();
    auto& file_info = param.get_file_info();
    auto& basic_info = param.get_basic_info();
    bool handled = false;

    auto pushUnique = [](std::vector<std::string>& v, const std::string& s) {
        if (std::find(v.begin(), v.end(), s) == v.end()) {
            v.push_back(s);
        }
    };

    if (key == "mctime" || key == "mc_time") {
        handled = true;
        // Support accumulation of multiple mctime lines
        auto times = InpParserStructures::parse_float_vector(value);
        for (float t : times) {
            mc_info.mc_time_list.push_back(t);
        }
    } else if (key == "energy_cutoff") {
        handled = true;
        // Legacy gcmc_opencl key: treat as a shared cutoff for fragment/protein.
        float cutoff = std::stof(value);
        if (cutoff > 0.0f) {
            energy_info.fragment_cutoff = cutoff;
            energy_info.fragment_cutoff_squared = cutoff * cutoff;
            energy_info.protein_cutoff = cutoff;
            energy_info.protein_cutoff_squared = cutoff * cutoff;
            // Mirror into the global cutoff so that energy calculations stay consistent.
            space_info.cutoff = cutoff;
            space_info.cutoff_explicit = true;
            space_info.cutoff_from_energy_cutoff = true;
        }
    } else if (key == "energy_cutoff_frag" || key == "energy_cutoff_fragment") {
        handled = true;
        float cutoff = std::stof(value);
        if (cutoff > 0.0f) {
            energy_info.fragment_cutoff = cutoff;
            energy_info.fragment_cutoff_squared = cutoff * cutoff;
            // Mirror fragment cutoff into the global cutoff so that energy calculations stay consistent
            space_info.cutoff = cutoff;
            space_info.cutoff_explicit = true;
            space_info.cutoff_from_energy_cutoff = true;
        }
    } else if (key == "energy_cutoff_prot" || key == "energy_cutoff_protein") {
        handled = true;
        float cutoff = std::stof(value);
        if (cutoff > 0.0f) {
            energy_info.protein_cutoff = cutoff;
            energy_info.protein_cutoff_squared = cutoff * cutoff;
            space_info.cutoff_explicit = true;
            space_info.cutoff_from_energy_cutoff = true;
        }
    } else if (key == "fragradius") {
        handled = true;
        // Raw value; normalized to internal nm in enhance_param.
        frag_info.radius_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "fragconf" || key == "fragconfs") {
        handled = true;
        // keep both for compatibility
        frag_info.conf_list = InpParserStructures::parse_int_vector(value);
        frag_info.fragconf_list = frag_info.conf_list;
    } else if (key == "num_conf_bias_trial" || key == "confbias_trials") {
        handled = true;
        bias_info.num_conf_bias_trials = static_cast<unsigned int>(std::stoi(value));
    } else if (key == "cavity_grid_dx" || key == "cavity_grid_spacing") {
        handled = true;
        // Map to grid_dx if given (fallback)
        // Raw value; normalized to internal nm in enhance_param.
        space_info.grid_spacing = std::stof(value);
    } else if (key == "cavity_grid_dx_frag") {
        handled = true;
        frag_info.cavity_grid_dx_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "probe_radius" || key == "cavity_probe_radius") {
        handled = true;
        // Map to sigma (approximate) if present
        // Raw value; normalized to internal nm in enhance_param.
        float r = std::stof(value);
        bias_info.sigma = r;
        bias_info.sigma_squared = r * r;
    } else if (key == "cavity_probe_radius_frag") {
        handled = true;
        frag_info.cavity_probe_radius_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "wdens") {
        handled = true;
        // Water density output control
        mc_info.wdens = std::stof(value);
    } else if (key == "eps" || key == "epsilon") {
        handled = true;
        // Dielectric constant
        frag_info.epsilon = std::stof(value);
    } else if (key == "target_numwaters" || key == "target_num_waters") {
        handled = true;
        // Target number of water molecules - store in both places
        frag_info.target_num_waters = std::stoi(value);
    } else if (key == "cavity_mask_frag") {
        handled = true;
        frag_info.cavity_mask_list = InpParserStructures::parse_int_vector(value);
    } else if (key == "gcmc_region") {
        handled = true;
        // GCMC insertion region (sphere/box specification)
        space_info.gcmc_region = value;
    } else if (key == "exclude_protein_volume") {
        handled = true;
        // Exclude protein volume from cavity bias
        space_info.exclude_protein_volume = (value == "yes" || value == "true" || value == "1");
    } else if (key == "use_vdw_radius_for_grid" || key == "use_vdw_radii_for_grid") {
        handled = true;
        // Use VDW radii for grid generation
        space_info.use_vdw_radius_for_grid = (value == "yes" || value == "true" || value == "1");
    } else if (key == "exclude_hydrogens_from_grid") {
        handled = true;
        // Exclude hydrogens from grid occupancy
        space_info.exclude_hydrogens_from_grid = (value == "yes" || value == "true" || value == "1");
    } else if (key == "use_switching") {
        handled = true;
        // Enable switching function
        mc_info.use_switching = (value == "yes" || value == "true" || value == "1");
    } else if (key == "switch_r_on" || key == "switch_ron") {
        handled = true;
        // Switching function r_on
        mc_info.switch_r_on = std::stof(value);  // Already in nm
    } else if (key == "switch_r_off" || key == "switch_roff") {
        handled = true;
        // Switching function r_off
        mc_info.switch_r_off = std::stof(value);  // Already in nm
    } else if (key == "switch_dist_frag" || key == "switch_dist_fragment") {
        handled = true;
        float dist = std::stof(value);
        energy_info.switch_dist_fragment = dist;
        energy_info.switch_dist_fragment_squared = dist * dist;
        energy_info.use_switching = true;
    } else if (key == "switch_dist_prot" || key == "switch_dist_protein") {
        handled = true;
        float dist = std::stof(value);
        energy_info.switch_dist_protein = dist;
        energy_info.switch_dist_protein_squared = dist * dist;
        energy_info.use_switching = true;
    } else if (key == "pairlist_freq") {
        handled = true;
        // Pairlist update frequency
        energy_info.pairlist_freq = static_cast<unsigned int>(std::stoi(value));
        // Pairlist rebuild scheduling is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "use_group_cutoff") {
        handled = true;
        // Use group-based cutoff instead of atom-based
        energy_info.use_group_cutoff = (value == "yes" || value == "true" || value == "1");
        // Group cutoff mode is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "pairlist_cutoff") {
        handled = true;
        // Pairlist cutoff distance for fragments
        energy_info.pairlist_cutoff = std::stof(value);  // Already in nm
        energy_info.pairlist_cutoff_squared = energy_info.pairlist_cutoff * energy_info.pairlist_cutoff;
        energy_info.pair_list_cutoff_fragment = energy_info.pairlist_cutoff;
        energy_info.pair_list_cutoff_fragment_squared = energy_info.pairlist_cutoff_squared;
        // Pairlist cutoff is not used by gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "pairlist_cutoff_protein") {
        handled = true;
        // Pairlist cutoff distance for protein
        float cutoff = std::stof(value);  // Already in nm
        energy_info.pair_list_cutoff_protein = cutoff;
        energy_info.pair_list_cutoff_protein_squared = cutoff * cutoff;
        // Protein pairlist cutoff is not used by gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "attempt_prob_ins") {
        handled = true;
        // Per-fragment insertion attempt probabilities
        mc_info.attempt_prob_ins = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_del") {
        handled = true;
        // Per-fragment deletion attempt probabilities
        mc_info.attempt_prob_del = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_trn") {
        handled = true;
        // Per-fragment translation attempt probabilities
        mc_info.attempt_prob_trn = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_rot") {
        handled = true;
        // Per-fragment rotation attempt probabilities
        mc_info.attempt_prob_rot = InpParserStructures::parse_float_vector(value);
    }

    (void)file_info;
    return handled;
}

} // namespace parameters
} // namespace io
} // namespace pygcmc

