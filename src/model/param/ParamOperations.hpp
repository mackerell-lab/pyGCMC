#pragma once

#ifndef PYGCMC_MODEL_PARAM_OPERATIONS_HPP
#define PYGCMC_MODEL_PARAM_OPERATIONS_HPP

#include "ParamStructures.hpp"
#include <stdexcept>
#include <cmath>

namespace pygcmc {
namespace model {
namespace param {

/**
 * @brief Static operations for parameter modification and updates
 */
class ParamOperations {
public:
    // === Update derived values ===
    static void updateBeta(MCParams& mc_params) {
        if (mc_params.temperature > 0.0f) {
            mc_params.beta = 1.0f / (mc_params.BOLTZMANN * mc_params.temperature);
        }
    }

    static void updateVolume(SpaceInfo& space_info) {
        space_info.volume = space_info.box_size[0] * space_info.box_size[1] * space_info.box_size[2];
    }

    static void updateEnergySquaredValues(EnergyInfo& energy_info) {
        energy_info.fragment_cutoff_squared = energy_info.fragment_cutoff * energy_info.fragment_cutoff;
        energy_info.protein_cutoff_squared = energy_info.protein_cutoff * energy_info.protein_cutoff;
        energy_info.pairlist_cutoff_squared = energy_info.pairlist_cutoff * energy_info.pairlist_cutoff;
        energy_info.switch_dist_fragment_squared = energy_info.switch_dist_fragment * energy_info.switch_dist_fragment;
        energy_info.switch_dist_protein_squared = energy_info.switch_dist_protein * energy_info.switch_dist_protein;
        energy_info.pair_list_cutoff_fragment_squared = energy_info.pair_list_cutoff_fragment * energy_info.pair_list_cutoff_fragment;
        energy_info.pair_list_cutoff_protein_squared = energy_info.pair_list_cutoff_protein * energy_info.pair_list_cutoff_protein;
    }

    static void updateFragmentSquaredValues(FragmentInfo& fragment_info) {
        fragment_info.init_cutoff_squared = fragment_info.init_cutoff * fragment_info.init_cutoff;
        fragment_info.gcmc_cutoff_squared = fragment_info.gcmc_cutoff * fragment_info.gcmc_cutoff;
    }

    static void updateBiasSquaredValues(BiasInfo& bias_info) {
        bias_info.sigma_squared = bias_info.sigma * bias_info.sigma;
    }

    // === Convenience setters ===
    static void setTemperature(MCParams& mc_params, float temperature) {
        if (temperature <= 0.0f) {
            throw std::invalid_argument("Temperature must be positive");
        }
        mc_params.temperature = temperature;
        updateBeta(mc_params);
    }

    static void setBoxSize(SpaceInfo& space_info, float x, float y, float z) {
        space_info.box_size[0] = x;
        space_info.box_size[1] = y;
        space_info.box_size[2] = z;
        updateVolume(space_info);
    }

    static void setMCSteps(MCParams& mc_params, int steps) {
        if (steps <= 0) {
            throw std::invalid_argument("MC steps must be positive");
        }
        mc_params.mc_steps = steps;
    }

    static void setCutoff(SpaceInfo& space_info, EnergyInfo& energy_info, float cutoff) {
        if (cutoff <= 0.0f) {
            throw std::invalid_argument("Cutoff must be positive");
        }
        space_info.cutoff = cutoff;
        energy_info.fragment_cutoff = cutoff;
        energy_info.protein_cutoff = cutoff;
        updateEnergySquaredValues(energy_info);
    }

    // === Clear operations ===
    static void clearBasicInfo(BasicInfo& basic_info) {
        basic_info = BasicInfo();
    }

    static void clearSpaceInfo(SpaceInfo& space_info) {
        space_info = SpaceInfo();
    }

    static void clearMCParams(MCParams& mc_params) {
        mc_params = MCParams();
    }

    static void clearEnergyInfo(EnergyInfo& energy_info) {
        energy_info = EnergyInfo();
    }

    static void clearFragmentInfo(FragmentInfo& fragment_info) {
        fragment_info = FragmentInfo();
    }

    static void clearBiasInfo(BiasInfo& bias_info) {
        bias_info = BiasInfo();
    }

    static void clearFileInfo(FileInfo& file_info) {
        file_info = FileInfo();
    }
};

} // namespace param
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_OPERATIONS_HPP