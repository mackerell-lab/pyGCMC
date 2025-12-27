#include "InpParserGCMC.hpp"
#include "../../model/param/ParamOperations.hpp"
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace io {
namespace parameters {

void InpParserGCMC::enhance_param_finalize(model::param::Param& param) {
    auto& mc_info = param.get_mc_info();
    auto& file_info = param.get_file_info();
    auto& energy_info = param.get_energy_info();
    auto& bias_info = param.get_bias_info();
    auto& space_info = param.get_space_info();
    auto& fragment_info = param.get_fragment_info();

    // Clamp insertion/deletion ratio into [0,1] and mirror to translation/rotation ratio
    mc_info.insertion_deletion_frac = std::max(0.0f, std::min(1.0f, mc_info.insertion_deletion_frac));
    mc_info.translation_rotation_frac = std::max(0.0f, 1.0f - mc_info.insertion_deletion_frac);

    // Ensure switching flags are synchronized between MC and energy configs
    if (mc_info.use_switching) {
        energy_info.use_switching = true;
        if (energy_info.switch_dist_fragment == 0.0f) {
            energy_info.switch_dist_fragment = mc_info.switch_r_on;
        }
        if (energy_info.switch_dist_protein == 0.0f) {
            energy_info.switch_dist_protein = mc_info.switch_r_off;
        }
    } else if (energy_info.use_switching) {
        mc_info.use_switching = true;
        if (mc_info.switch_r_on <= 0.0f) {
            mc_info.switch_r_on = energy_info.switch_dist_fragment;
        }
        if (mc_info.switch_r_off <= 0.0f) {
            mc_info.switch_r_off = energy_info.switch_dist_protein > 0.0f ?
                energy_info.switch_dist_protein : energy_info.switch_dist_fragment;
        }
    }

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

    // If no per-fragment move probabilities were provided, derive them from insdel_frac
    if (mc_info.attempt_prob_ins.empty() &&
        mc_info.attempt_prob_del.empty() &&
        mc_info.attempt_prob_trn.empty() &&
        mc_info.attempt_prob_rot.empty()) {

        size_t fragment_count = file_info.fragment_names.empty() ? 1 : file_info.fragment_names.size();
        const float ins_prob = mc_info.insertion_deletion_frac * 0.5f;
        const float del_prob = mc_info.insertion_deletion_frac * 0.5f;
        const float trn_prob = mc_info.translation_rotation_frac * 0.5f;
        const float rot_prob = mc_info.translation_rotation_frac * 0.5f;

        mc_info.attempt_prob_ins.assign(fragment_count, ins_prob);
        mc_info.attempt_prob_del.assign(fragment_count, del_prob);
        mc_info.attempt_prob_trn.assign(fragment_count, trn_prob);
        mc_info.attempt_prob_rot.assign(fragment_count, rot_prob);
    }

    const size_t fragment_count = file_info.fragment_names.empty() ? 1 : file_info.fragment_names.size();
    auto ensureFloatList = [&](std::vector<float>& vec, float fallback) {
        if (vec.empty()) {
            vec.assign(fragment_count, fallback);
        } else if (vec.size() < fragment_count) {
            vec.resize(fragment_count, fallback);
        }
    };
    auto ensureIntList = [&](std::vector<int>& vec, int fallback) {
        if (vec.empty()) {
            vec.assign(fragment_count, fallback);
        } else if (vec.size() < fragment_count) {
            vec.resize(fragment_count, fallback);
        }
    };

    const float defaultGrid = space_info.grid_spacing > 0.0f ? space_info.grid_spacing : 0.2f;
    const float defaultProbe = bias_info.sigma > 0.0f ? bias_info.sigma : 0.14f;
    ensureFloatList(fragment_info.cavity_grid_dx_list, defaultGrid);
    ensureFloatList(fragment_info.cavity_probe_radius_list, defaultProbe);
    ensureIntList(fragment_info.cavity_mask_list, -1);

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

    // Update derived squared quantities (safe to repeat; units conversion is idempotent)
    model::param::ParamOperations::updateEnergySquaredValues(energy_info);
    model::param::ParamOperations::updateFragmentSquaredValues(fragment_info);
    model::param::ParamOperations::updateBiasSquaredValues(bias_info);
}

} // namespace parameters
} // namespace io
} // namespace pygcmc

