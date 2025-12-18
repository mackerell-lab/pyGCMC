#pragma once

#ifndef PYGCMC_MODEL_PARAM_STRUCTURES_HPP
#define PYGCMC_MODEL_PARAM_STRUCTURES_HPP

#include <vector>
#include <string>
#include <array>

namespace pygcmc {
namespace model {
namespace param {

/**
 * @brief Basic information structure for GCMC simulation
 */
struct BasicInfo {
    std::string version = "gcmc_v2.0";
    int verbosity = 0;
    int debug = 0;
    bool print_logfile = false;
    std::string param_file;
    std::string log_file;
    // INP unit system for compatibility with legacy gcmc_gpu style inputs.
    // Supported values (case-insensitive):
    // - "auto" (default) => assume gcmc_gpu/opencl style (Å + kcal/mol)
    // - "nm" (alias: "openmm") => native/internal nm + kJ/mol decks
    // - "gcmc_gpu"/"charmm"/"a"/"angstrom" => Å + kcal/mol decks
    std::string inp_units = "auto";
    // Whether the user explicitly specified inp_units/units in the INP.
    // This is used to avoid overriding user intent when applying legacy version heuristics.
    bool inp_units_explicit = false;
    // Internal flag to make enhance_param() idempotent for unit conversion.
    bool inp_units_converted = false;
    // Diagnostics: capture INP keys seen/handled/unknown during parsing.
    // These are used to avoid silently ignoring legacy keys during gcmc_gpu/opencl compatibility work.
    std::vector<std::string> inp_keys_seen;
    std::vector<std::string> inp_keys_handled;
    std::vector<std::string> inp_keys_unknown;
    unsigned int random_seed = 0;
    int num_threads = 1;
    bool is_box = false;
    bool init_cycle = false;
    bool conserve_fragments = false;
};

/**
 * @brief System space structure for simulation geometry
 */
struct SpaceInfo {
    float grid_spacing = 1.0;
    std::array<float, 3> gc_center = {0.0, 0.0, 0.0};
    std::array<float, 3> sys_center = {0.0, 0.0, 0.0};
    std::array<float, 3> crystal_dim = {0.0, 0.0, 0.0};
    std::array<float, 3> box_size = {0.0, 0.0, 0.0};

    float volume = 0.0;
    float target_volume = 0.0;
    float sys_box_volume = 0.0;
    float gcmc_volume = 0.0;
    float protein_volume = 0.0;
    
    bool use_vdw_radius_for_grid = false;
    bool exclude_hydrogens_from_grid = false;
    bool exclude_protein_volume = false;

    std::string gcmc_region = "";  // GCMC insertion region specification

    float tmp_prob = 0.0;
    float cutoff = 12.0;
    bool cutoff_explicit = false;
    // True when the cutoff was provided via legacy gcmc_opencl keys (energy_cutoff*),
    // which are interpreted as Å in inp_units:auto mode.
    bool cutoff_from_energy_cutoff = false;
};

/**
 * @brief Monte Carlo simulation parameters
 */
struct MCParams {
    int mc_steps = 1;
    int moves_per_step = 1;  // Number of GCMC moves per MC step
    int current_step = 0;
    int print_freq = 1;
    int save_freq = 10000;  // nsave: trajectory save frequency

    float temperature = 300.0;
    float beta = 1.0;

    float insertion_deletion_frac = 0.5;
    float translation_rotation_frac = 0.5;

    float max_translation_dist = 1.0;
    float max_rotation_angle = 30.0;

    std::vector<std::string> operation_types = {"Ins", "Del", "Trn", "Rot"};
    std::vector<float> mc_time_list;
    std::vector<float> mc_time_cumulative;

    // Per-fragment move probabilities (legacy compatibility)
    std::vector<float> attempt_prob_ins;   // Insert attempt probabilities per fragment
    std::vector<float> attempt_prob_del;   // Delete attempt probabilities per fragment
    std::vector<float> attempt_prob_trn;   // Translate attempt probabilities per fragment
    std::vector<float> attempt_prob_rot;   // Rotate attempt probabilities per fragment

    // Analysis and control
    float wdens = 0.0f;  // Water density value or output frequency

    // Switching function parameters
    bool use_switching = false;
    float switch_r_on = 10.0;  // Switching function start distance
    float switch_r_off = 12.0; // Switching function cutoff distance

    std::vector<float> fragment_prob;
    std::vector<float> water_prob;
    std::vector<float> atom_prob;
    std::vector<float> test_prob;

    int rotate_dih_status = 0;

    // Physical constants (use kJ/mol/K to match platform energy units)
    float BOLTZMANN = 8.314e-3f;
    float KCAL_TO_KJ = 4.184f;
};

/**
 * @brief Energy calculation parameters
 */
struct EnergyInfo {
    bool use_group_cutoff = true;
    float fragment_cutoff = 10.0;
    float protein_cutoff = 10.0;
    float fragment_cutoff_squared = 100.0;
    float protein_cutoff_squared = 100.0;

    float pairlist_cutoff = 0.0;
    float pairlist_cutoff_squared = 0.0;
    unsigned int pairlist_freq = 1000;

    bool use_switching = false;
    float switch_dist_fragment = 0.0;
    float switch_dist_protein = 0.0;
    float switch_dist_fragment_squared = 0.0;
    float switch_dist_protein_squared = 0.0;

    float energy_sw_ref = 1.0;
    float energy_sw_scale = 1.0;

    bool test_sw_filters = false;
    bool apply_sw_filters = false;
    bool test_energy = false;

    float pair_list_cutoff_fragment = 0.0;
    float pair_list_cutoff_protein = 0.0;
    float pair_list_cutoff_fragment_squared = 0.0;
    float pair_list_cutoff_protein_squared = 0.0;
};

/**
 * @brief Fragment parameters for GCMC simulation
 */
struct FragmentInfo {
    float water_density = 55.0;
    float epsilon = 1.0;
    int num_waters = 0;
    int target_num_waters = 0;
    int water_index = 0;

    float excess_threshold = 1.0;

    bool use_number_water_nbar = false;
    bool use_const_water_nbar = false;
    int const_water_nbar = 0;

    float init_cutoff = 0.0;
    float init_cutoff_squared = 0.0;
    bool use_gcmc_cutoff = false;
    float gcmc_cutoff = 0.0;
    float gcmc_cutoff_squared = 0.0;

    std::vector<int> remove_init;
    std::vector<int> remove_excess;

    std::vector<int> confs_list;
    std::vector<int> cavity_index_list;
    std::vector<float> cavity_list;
    std::vector<float> cavity_grid_dx_list;
    std::vector<float> cavity_probe_radius_list;
    std::vector<int> cavity_mask_list;

    std::vector<float> conc_list;
    std::vector<float> muex_list;
    std::vector<float> radius_list;

    std::vector<int> conf_list;
    int flag_remove_init = 0;
    int flag_remove_excess = 0;
    int total_protitp_size = 0;
    std::vector<int> fragconf_list;
};

/**
 * @brief Biased sampling parameters
 */
struct BiasInfo {
    bool use_cavity_bias = false;
    float sigma = 2.4;
    float sigma_squared = 5.76;

    bool use_conf_bias = false;
    unsigned int num_conf_bias_trials = 10;
};

/**
 * @brief File path parameters
 */
struct FileInfo {
    std::string topology_file;
    std::string input_pdb_file;
    std::string output_pdb_file;
    std::string output_top_file;

    std::string atomtype_file;
    std::string monomer_dir;
    std::string conc_norm = "water";
    std::string conc_region = "total";
    std::vector<std::string> par_files;

    std::vector<std::string> protein_top_files;
    std::vector<std::string> fragment_top_files;
    std::vector<std::string> fragment_names;
    std::vector<std::string> fragment_mqtr_files;
    std::string tmp_frag_name;

    bool generate_maps = false;
    std::string map_prefix = "gc_maps";
};

} // namespace param
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_STRUCTURES_HPP
