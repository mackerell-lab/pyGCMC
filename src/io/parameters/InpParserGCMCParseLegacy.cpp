#include "InpParserGCMC.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

namespace pygcmc {
namespace io {
namespace parameters {

using ::pygcmc::io::parameters::InpParserStructures;

bool InpParserGCMC::parse_line_ext_legacy(const std::string& key, const std::string& value, model::param::Param& param) {
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

    if (key == "mc_move_prob") {
        handled = true;
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
        handled = true;
        // Fixed target number of water molecules
        frag_info.use_const_water_nbar = true;
        frag_info.const_water_nbar = static_cast<int>(std::stoi(value));
    } else if (key == "number_water_nbar") {
        handled = true;
        // Toggle number-based nbar (use current water count as target)
        frag_info.use_number_water_nbar = (value == "yes" || value == "true" || value == "1");
        if (frag_info.use_number_water_nbar) {
            frag_info.use_const_water_nbar = false;
        }
    } else if (key == "volume_water_nbar") {
        handled = true;
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
    } else if (key == "max_translation" || key == "max_translation_dist") {
        handled = true;
        // Raw value; normalized to internal nm in enhance_param.
        mc_info.max_translation_dist = std::stof(value);
    } else if (key == "max_rotation" || key == "max_rotation_angle") {
        handled = true;
        // Rotation angles are treated as degrees in legacy INP decks (converted to radians where needed).
        mc_info.max_rotation_angle = std::stof(value);
    } else if (key == "insdel_frac" || key == "insdel_fraction") {
        handled = true;
        float frac = std::stof(value);
        mc_info.insertion_deletion_frac = std::max(0.0f, std::min(1.0f, frac));
        mc_info.translation_rotation_frac = 1.0f - mc_info.insertion_deletion_frac;
    } else if (key == "attempt_prob_frag") {
        handled = true;
        mc_info.fragment_prob = InpParserStructures::parse_float_vector(value);
        // Legacy compatibility key; currently not used by gcmc_cpu (move scheduling uses attempt_prob_*).
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "attempt_prob_water") {
        handled = true;
        mc_info.water_prob = InpParserStructures::parse_float_vector(value);
        // Legacy compatibility key; currently not used by gcmc_cpu.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "attempt_prob_atom") {
        handled = true;
        mc_info.atom_prob = InpParserStructures::parse_float_vector(value);
        // Legacy compatibility key; currently not used by gcmc_cpu.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "test_energy") {
        handled = true;
        energy_info.test_energy = (value == "yes" || value == "true" || value == "1");
        // Energy diagnostic mode is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "test_sw_filters" || key == "test_SW_filters") {
        handled = true;
        energy_info.test_sw_filters = (value == "yes" || value == "true" || value == "1");
        // Switching filter diagnostics are not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "apply_sw_filters" || key == "apply_SW_filters") {
        handled = true;
        energy_info.apply_sw_filters = (value == "yes" || value == "true" || value == "1");
        // Switching filter application is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "sw_reference" || key == "SW_reference") {
        handled = true;
        // Legacy inputs provide kcal/mol – convert to kJ/mol for internal use
        energy_info.energy_sw_ref = std::stof(value) * 4.184f;
        // Switching filter reference is not used by gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "sw_scale" || key == "SW_scale") {
        handled = true;
        energy_info.energy_sw_scale = std::stof(value) * 4.184f;
        // Switching filter scale is not used by gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "rotate_dihedral" || key == "rotate_dih_status") {
        handled = true;
        // Legacy switch: enable/disable dihedral rotation
        if (value == "yes" || value == "true" || value == "1") {
            mc_info.rotate_dih_status = 1;
        } else if (value == "no" || value == "false" || value == "0") {
            mc_info.rotate_dih_status = 0;
        } else {
            // Allow numeric value passthrough
            try {
                mc_info.rotate_dih_status = std::stoi(value);
            } catch (const std::exception&) {
                mc_info.rotate_dih_status = 0;
            }
        }
        // Dihedral rotation moves are not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "remove_init") {
        handled = true;
        frag_info.remove_init = InpParserStructures::parse_int_vector(value);
        frag_info.flag_remove_init = 1;
        // Initial-fragment removal logic is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "remove_excess") {
        handled = true;
        frag_info.remove_excess = InpParserStructures::parse_int_vector(value);
        frag_info.flag_remove_excess = 1;
        // Excess-fragment removal logic is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "initial_fragments_cutoff") {
        handled = true;
        const float cutoff = std::stof(value);
        frag_info.init_cutoff = cutoff;
        frag_info.init_cutoff_squared = cutoff * cutoff;
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "excess_fragments_threshold") {
        handled = true;
        frag_info.excess_threshold = std::stof(value);
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "gcmc_cutoff") {
        handled = true;
        frag_info.gcmc_cutoff = std::stof(value);
        frag_info.gcmc_cutoff_squared = frag_info.gcmc_cutoff * frag_info.gcmc_cutoff;
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "use_gcmc_cutoff") {
        handled = true;
        frag_info.use_gcmc_cutoff = (value == "yes" || value == "true" || value == "1");
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "target_volume") {
        handled = true;
        space_info.target_volume = std::stof(value);
        // Target volume is not used by gcmc_cpu yet (stored for compatibility only).
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "use_const_water_nbar") {
        handled = true;
        // gcmc_gpu compatibility: allow either yes/no or an integer value
        try {
            const int n = std::stoi(value);
            if (n > 0) {
                frag_info.use_const_water_nbar = true;
                frag_info.const_water_nbar = n;
            } else {
                frag_info.use_const_water_nbar = false;
            }
        } catch (const std::exception&) {
            frag_info.use_const_water_nbar = (value == "yes" || value == "true" || value == "1");
        }
    } else if (key == "use_number_water_nbar") {
        handled = true;
        frag_info.use_number_water_nbar = (value == "yes" || value == "true" || value == "1");
        if (frag_info.use_number_water_nbar) {
            frag_info.use_const_water_nbar = false;
        }
    } else if (key == "fragmqtr") {
        handled = true;
        // Legacy: per-fragment MQTR file(s)
        file_info.fragment_mqtr_files.push_back(value);
        // MQTR functionality is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "inp_units" || key == "units") {
        handled = true;
        // Override unit system ("auto", "nm", "gcmc_gpu"/"angstrom"/"a")
        std::string v = value;
        std::transform(v.begin(), v.end(), v.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        basic_info.inp_units = v;
        basic_info.inp_units_explicit = true;
    }

    (void)bias_info;
    return handled;
}

} // namespace parameters
} // namespace io
} // namespace pygcmc

