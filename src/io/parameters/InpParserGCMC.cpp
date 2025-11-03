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
        frag_info.radius_list = InpParserStructures::parse_float_vector(value);  // Already in nm
    } else if (key == "fragconf" || key == "fragconfs") {
        // keep both for compatibility
        frag_info.conf_list = InpParserStructures::parse_int_vector(value);
        frag_info.fragconf_list = frag_info.conf_list;
    } else if (key == "num_conf_bias_trial" || key == "confbias_trials") {
        bias_info.num_conf_bias_trials = static_cast<unsigned int>(std::stoi(value));
    } else if (key == "cavity_grid_dx") {
        // Map to grid_dx if given (fallback)
        space_info.grid_spacing = std::stof(value);  // Already in nm
    } else if (key == "probe_radius") {
        // Map to sigma (approximate) if present
        float r = std::stof(value);  // Already in nm
        bias_info.sigma = r;
        bias_info.sigma_squared = r * r;
    } else if (key == "wdens") {
        // Water density output control
        mc_info.wdens = std::stof(value);
    } else if (key == "eps" || key == "epsilon") {
        // Dielectric constant
        frag_info.epsilon = std::stof(value);
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
        mc_info.switch_r_on = std::stof(value);  // Already in nm
    } else if (key == "switch_r_off" || key == "switch_roff") {
        // Switching function r_off
        mc_info.switch_r_off = std::stof(value);  // Already in nm
    } else if (key == "pairlist_freq") {
        // Pairlist update frequency
        energy_info.pairlist_freq = static_cast<unsigned int>(std::stoi(value));
    } else if (key == "use_group_cutoff") {
        // Use group-based cutoff instead of atom-based
        energy_info.use_group_cutoff = (value == "yes" || value == "true" || value == "1");
    } else if (key == "pairlist_cutoff") {
        // Pairlist cutoff distance for fragments
        energy_info.pairlist_cutoff = std::stof(value);  // Already in nm
        energy_info.pairlist_cutoff_squared = energy_info.pairlist_cutoff * energy_info.pairlist_cutoff;
        energy_info.pair_list_cutoff_fragment = energy_info.pairlist_cutoff;
        energy_info.pair_list_cutoff_fragment_squared = energy_info.pairlist_cutoff_squared;
    } else if (key == "pairlist_cutoff_protein") {
        // Pairlist cutoff distance for protein
        float cutoff = std::stof(value);  // Already in nm
        energy_info.pair_list_cutoff_protein = cutoff;
        energy_info.pair_list_cutoff_protein_squared = cutoff * cutoff;
    } else if (key == "attempt_prob_ins") {
        // Per-fragment insertion attempt probabilities
        mc_info.attempt_prob_ins = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_del") {
        // Per-fragment deletion attempt probabilities
        mc_info.attempt_prob_del = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_trn") {
        // Per-fragment translation attempt probabilities
        mc_info.attempt_prob_trn = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_rot") {
        // Per-fragment rotation attempt probabilities
        mc_info.attempt_prob_rot = InpParserStructures::parse_float_vector(value);
    } else if (key == "mc_move_prob") {
        // Legacy format: four weights [insert, delete, translate, rotate]
        // Store in attempt_prob_* temporarily, will be broadcasted in enhance_param
        auto probs = InpParserStructures::parse_float_vector(value);

        if (probs.size() >= 4) {
            // Clear any existing values and store the four probabilities
            // These will be broadcasted to all fragments in enhance_param()
            mc_info.attempt_prob_ins = {probs[0]};
            mc_info.attempt_prob_del = {probs[1]};
            mc_info.attempt_prob_trn = {probs[2]};
            mc_info.attempt_prob_rot = {probs[3]};

            // Set a flag to indicate mc_move_prob was used (store in first element as negative to mark)
            // This is a hack to avoid adding new fields to MCInfo
            // We'll check this in enhance_param and broadcast to all fragments
            std::cout << "[INP] Parsed mc_move_prob: "
                      << probs[0] << " (ins), "
                      << probs[1] << " (del), "
                      << probs[2] << " (trn), "
                      << probs[3] << " (rot)"
                      << " - will be applied to all fragments" << std::endl;
        } else {
            std::cerr << "[WARNING] mc_move_prob requires 4 values, got "
                      << probs.size() << std::endl;
        }
    } else if (key == "const_water_nbar") {
        // Fixed target number of water molecules
        frag_info.use_const_water_nbar = true;
        frag_info.const_water_nbar = static_cast<int>(std::stoi(value));
    } else if (key == "number_water_nbar") {
        // Toggle number-based nbar (use current water count as target)
        frag_info.use_number_water_nbar = (value == "yes" || value == "true" || value == "1");
        if (frag_info.use_number_water_nbar) {
            frag_info.use_const_water_nbar = false;
        }
    } else if (key == "volume_water_nbar") {
        // Volume-based nbar via target concentration (M)
        // Use this to override water concentration for volume projection mode
        // (default mode already uses concentration, but we record explicitly)
        const float concM = std::stof(value);
        if (!frag_info.conc_list.empty()) {
            // Override first (water) entry if present
            frag_info.conc_list[0] = concM;
        } else {
            frag_info.conc_list.push_back(concM);
        }
    }
}

void InpParserGCMC::enhance_param(model::param::Param& param) {
    auto& mc_info = param.get_mc_info();
    auto& file_info = param.get_file_info();

    // Broadcast mc_move_prob to all fragments if it was used
    // If attempt_prob_* vectors have size 1, it means mc_move_prob was set
    // We need to broadcast to all fragments
    if (mc_info.attempt_prob_ins.size() == 1 &&
        mc_info.attempt_prob_del.size() == 1 &&
        mc_info.attempt_prob_trn.size() == 1 &&
        mc_info.attempt_prob_rot.size() == 1) {

        size_t fragment_count = file_info.fragment_names.size();
        if (fragment_count > 1) {
            // Broadcast the single values to all fragments
            float ins_val = mc_info.attempt_prob_ins[0];
            float del_val = mc_info.attempt_prob_del[0];
            float trn_val = mc_info.attempt_prob_trn[0];
            float rot_val = mc_info.attempt_prob_rot[0];

            mc_info.attempt_prob_ins.assign(fragment_count, ins_val);
            mc_info.attempt_prob_del.assign(fragment_count, del_val);
            mc_info.attempt_prob_trn.assign(fragment_count, trn_val);
            mc_info.attempt_prob_rot.assign(fragment_count, rot_val);

            std::cout << "[INP] Broadcasted mc_move_prob to " << fragment_count << " fragments" << std::endl;
        }
    }

    // Ensure mc_time_cumulative if mc_time_list provided
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