#pragma once

#ifndef PYGCMC_MODEL_PARAM_MAIN_HPP
#define PYGCMC_MODEL_PARAM_MAIN_HPP

#include "../common/ModelInterface.hpp"
#include "../common/ModelConstants.hpp"
#include <vector>
#include <string>
#include <array>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace pygcmc {
namespace model {
namespace param {

/**
 * @brief Basic information structure for GCMC simulation
 */
struct BasicInfo {
    std::string version = "gcmc_v2.0";
    int verbosity = 0;                        ///< Output verbosity level (0-2)
    int debug = 0;                           ///< Debug level (0-2)
    bool print_logfile = false;              ///< Whether to output log file
    std::string param_file;                  ///< Parameter file path
    std::string log_file;                    ///< Log file path
    unsigned int random_seed = 0;            ///< Random number seed
    int num_threads = 1;                     ///< Number of parallel threads
    bool is_box = false;                     ///< Box size flag
    bool init_cycle = false;                 ///< Whether to start from initial configuration
    bool conserve_fragments = false;         ///< Whether to fix solute quantity
};

/**
 * @brief System space structure for simulation geometry
 */
struct SpaceInfo {
    float grid_spacing = 1.0;                ///< Grid spacing for cavity detection (Å)
    std::array<float, 3> gc_center = {0.0, 0.0, 0.0};    ///< GCMC region center
    std::array<float, 3> sys_center = {0.0, 0.0, 0.0};   ///< System center coordinates
    std::array<float, 3> crystal_dim = {0.0, 0.0, 0.0};  ///< Crystal cell dimensions
    std::array<float, 3> box_size = {0.0, 0.0, 0.0};     ///< Simulation box dimensions (Å)

    float volume = 0.0;                      ///< Total volume
    float target_volume = 0.0;               ///< Target volume
    float sys_box_volume = 0.0;              ///< System box volume
    float gcmc_volume = 0.0;                 ///< GCMC region volume
    float protein_volume = 0.0;              ///< Protein volume
    
    bool use_vdw_radius_for_grid = false;    ///< Use van der Waals radius for grid
    bool exclude_hydrogens_from_grid = false; ///< Exclude hydrogens from grid
    bool exclude_protein_volume = false;      ///< Exclude protein volume
    
    float tmp_prob = 0.0;                    ///< Temporary probability variable
    float cutoff = 12.0;                     ///< Non-bonded interaction cutoff (Å)

    /**
     * @brief Calculate volume from box dimensions
     */
    void updateVolume() {
        volume = box_size[0] * box_size[1] * box_size[2];
    }

    /**
     * @brief Set box dimensions and update volume
     */
    void setBoxSize(float x, float y, float z) {
        box_size[0] = x;
        box_size[1] = y;
        box_size[2] = z;
        updateVolume();
    }
};

/**
 * @brief Monte Carlo simulation parameters
 */
struct MCParams {
    int mc_steps = 1;                        ///< Total MC steps
    int current_step = 0;                    ///< Current MC step
    int print_freq = 1;                      ///< Output frequency

    float temperature = 300.0;               ///< Simulation temperature (K)
    float beta = 1.0;                        ///< β = 1/(kB*T)

    float insertion_deletion_frac = 0.5;     ///< Insertion/deletion operation ratio
    float translation_rotation_frac = 0.5;    ///< Translation/rotation operation ratio

    float max_translation_dist = 1.0;        ///< Maximum translation distance (Å)
    float max_rotation_angle = 30.0;         ///< Maximum rotation angle (degrees)

    std::vector<std::string> operation_types = {"Ins", "Del", "Trn", "Rot"};  ///< MC operation types
    std::vector<float> mc_time_list;         ///< Time weights for move types
    std::vector<float> mc_time_cumulative;   ///< Cumulative time for selection

    std::vector<float> fragment_prob;        ///< Fragment operation probabilities
    std::vector<float> water_prob;           ///< Water molecule operation probabilities
    std::vector<float> atom_prob;            ///< Atom operation probabilities
    std::vector<float> test_prob;            ///< Test operation probabilities

    int rotate_dih_status = 0;               ///< Dihedral rotation status

    // Physical constants
    float BOLTZMANN = 0.001987f;             ///< Boltzmann constant (kcal/mol/K)
    float KCAL_TO_KJ = 4.184f;               ///< Energy unit conversion factor

    /**
     * @brief Update beta from temperature
     */
    void updateBeta() {
        if (temperature > 0.0f) {
            beta = 1.0f / (BOLTZMANN * temperature);
        }
    }

    /**
     * @brief Set temperature and update beta
     */
    void setTemperature(float temp) {
        if (temp <= 0.0f) {
            throw std::invalid_argument("Temperature must be positive");
        }
        temperature = temp;
        updateBeta();
    }
};

/**
 * @brief Energy calculation parameters
 */
struct EnergyInfo {
    bool use_group_cutoff = true;            ///< Use group cutoff optimization
    float fragment_cutoff = 10.0;            ///< Fragment energy cutoff (Å)
    float protein_cutoff = 10.0;             ///< Protein energy cutoff (Å)
    float fragment_cutoff_squared = 100.0;   ///< Square of fragment cutoff
    float protein_cutoff_squared = 100.0;    ///< Square of protein cutoff

    float pairlist_cutoff = 0.0;            ///< Pairlist cutoff (Å)
    float pairlist_cutoff_squared = 0.0;     ///< Square of pairlist cutoff
    unsigned int pairlist_freq = 1000;       ///< Pairlist update frequency

    bool use_switching = false;              ///< Use switching function
    float switch_dist_fragment = 0.0;        ///< Fragment switching distance (Å)
    float switch_dist_protein = 0.0;         ///< Protein switching distance (Å)
    float switch_dist_fragment_squared = 0.0; ///< Square of fragment switching distance
    float switch_dist_protein_squared = 0.0;  ///< Square of protein switching distance

    float energy_sw_ref = 1.0;               ///< SW reference energy
    float energy_sw_scale = 1.0;             ///< SW energy scale

    bool test_sw_filters = false;            ///< Test SW filters
    bool apply_sw_filters = false;           ///< Apply SW filters
    bool test_energy = false;                ///< Test energy

    float pair_list_cutoff_fragment = 0.0;   ///< Fragment pairlist cutoff (Å)
    float pair_list_cutoff_protein = 0.0;    ///< Protein pairlist cutoff (Å)
    float pair_list_cutoff_fragment_squared = 0.0;  ///< Square of fragment pairlist cutoff
    float pair_list_cutoff_protein_squared = 0.0;   ///< Square of protein pairlist cutoff

    /**
     * @brief Update squared cutoff values
     */
    void updateSquaredValues() {
        fragment_cutoff_squared = fragment_cutoff * fragment_cutoff;
        protein_cutoff_squared = protein_cutoff * protein_cutoff;
        pairlist_cutoff_squared = pairlist_cutoff * pairlist_cutoff;
        switch_dist_fragment_squared = switch_dist_fragment * switch_dist_fragment;
        switch_dist_protein_squared = switch_dist_protein * switch_dist_protein;
        pair_list_cutoff_fragment_squared = pair_list_cutoff_fragment * pair_list_cutoff_fragment;
        pair_list_cutoff_protein_squared = pair_list_cutoff_protein * pair_list_cutoff_protein;
    }
};

/**
 * @brief Fragment parameters for GCMC simulation
 */
struct FragmentInfo {
    float water_density = 55.0;              ///< Target water density (M)
    float epsilon = 1.0;                     ///< Epsilon parameter
    int num_waters = 0;                      ///< Current number of water molecules
    int target_num_waters = 0;               ///< Target number of water molecules
    int water_index = 0;                     ///< Index of water fragment

    float excess_threshold = 1.0;            ///< Excess threshold (L parameter)

    bool use_number_water_nbar = true;       ///< Use water molecule number mean
    bool use_const_water_nbar = false;       ///< Use constant water number mean
    int const_water_nbar = 0;                ///< Constant water number mean

    float init_cutoff = 0.0;                 ///< Initialization cutoff distance (Å)
    float init_cutoff_squared = 0.0;         ///< Square of initialization cutoff
    bool use_gcmc_cutoff = false;            ///< Use GCMC cutoff
    float gcmc_cutoff = 0.0;                 ///< GCMC region cutoff distance (Å)
    float gcmc_cutoff_squared = 0.0;         ///< Square of GCMC cutoff

    std::vector<int> remove_init;            ///< Initial removal list
    std::vector<int> remove_excess;          ///< Excess removal list

    std::vector<int> confs_list;             ///< Available configuration list
    std::vector<int> cavity_index_list;      ///< Cavity index list
    std::vector<float> cavity_list;          ///< Cavity list

    std::vector<float> conc_list;            ///< Target concentration list (M)
    std::vector<float> muex_list;            ///< Excess chemical potential list (kcal/mol)
    std::vector<float> radius_list;          ///< Molecule radius list (Å)

    std::vector<int> conf_list;              ///< Configuration list
    int flag_remove_init = 0;                ///< Initial removal flag
    int flag_remove_excess = 0;              ///< Excess removal flag
    int total_protitp_size = 0;              ///< Total number of protein topology files
    std::vector<int> fragconf_list;          ///< Fragment configuration list

    /**
     * @brief Update squared cutoff values
     */
    void updateSquaredValues() {
        init_cutoff_squared = init_cutoff * init_cutoff;
        gcmc_cutoff_squared = gcmc_cutoff * gcmc_cutoff;
    }
};

/**
 * @brief Biased sampling parameters
 */
struct BiasInfo {
    bool use_cavity_bias = false;            ///< Use cavity-biased sampling
    float sigma = 2.4;                       ///< Cavity size parameter (Å)
    float sigma_squared = 5.76;              ///< Square of sigma

    bool use_conf_bias = false;              ///< Use configuration-biased sampling
    unsigned int num_conf_bias_trials = 10;   ///< Configuration bias trial count

    /**
     * @brief Update squared sigma value
     */
    void updateSquaredValues() {
        sigma_squared = sigma * sigma;
    }
};

/**
 * @brief File path parameters
 */
struct FileInfo {
    std::string topology_file;               ///< System topology file
    std::string input_pdb_file;              ///< Initial configuration file
    std::string output_pdb_file;             ///< Trajectory output file
    std::string output_top_file;             ///< Output topology file

    std::string atomtype_file;               ///< Atom type definition file
    std::string monomer_dir;                 ///< Monomer directory
    std::string conc_norm = "water";         ///< Concentration normalization method
    std::string conc_region = "total";       ///< Concentration calculation region
    std::vector<std::string> par_files;      ///< Force field parameter files

    std::vector<std::string> protein_top_files;  ///< Protein topology files
    std::vector<std::string> fragment_top_files; ///< Fragment topology files
    std::vector<std::string> fragment_names;     ///< Fragment name list
    std::vector<std::string> fragment_mqtr_files;///< Fragment MQTR files
    std::string tmp_frag_name;               ///< Temporary fragment name

    bool generate_maps = false;              ///< Generate mapping files
    std::string map_prefix = "gc_maps";      ///< Mapping file name prefix

    /**
     * @brief Get fragment coordinate filename
     */
    std::string get_fragment_coordinate_filename(int frag_index) const {
        if (frag_index < 0 || frag_index >= static_cast<int>(fragment_names.size())) {
            throw std::out_of_range("Fragment index out of range");
        }
        return monomer_dir + "/" + fragment_names[frag_index] + ".pdb";
    }

    /**
     * @brief Check if all required files are specified
     */
    bool has_required_files() const {
        return !topology_file.empty() && !input_pdb_file.empty() && !output_pdb_file.empty();
    }
};

/**
 * @brief Complete GCMC parameter class with full functionality and backward compatibility
 * This class maintains the same API as the original param.hpp
 */
class Param : public common::IValidatable {
public:
    Param() = default;
    ~Param() = default;

    // IValidatable interface
    bool is_valid() const override {
        // Check basic validity
        if (mc_info_.temperature <= 0.0f) return false;
        if (mc_info_.mc_steps <= 0) return false;
        if (!file_info_.has_required_files()) return false;
        
        // Check energy parameters
        if (energy_info_.fragment_cutoff <= 0.0f || energy_info_.protein_cutoff <= 0.0f) return false;
        
        return true;
    }

    // Getters (const)
    const BasicInfo& get_basic_info() const { return basic_info_; }
    const SpaceInfo& get_space_info() const { return space_info_; }
    const MCParams& get_mc_info() const { return mc_info_; }
    const EnergyInfo& get_energy_info() const { return energy_info_; }
    const FragmentInfo& get_fragment_info() const { return fragment_info_; }
    const BiasInfo& get_bias_info() const { return bias_info_; }
    const FileInfo& get_file_info() const { return file_info_; }

    // Non-const getters for modification
    BasicInfo& get_basic_info() { return basic_info_; }
    SpaceInfo& get_space_info() { return space_info_; }
    MCParams& get_mc_info() { return mc_info_; }
    EnergyInfo& get_energy_info() { return energy_info_; }
    FragmentInfo& get_fragment_info() { return fragment_info_; }
    BiasInfo& get_bias_info() { return bias_info_; }
    FileInfo& get_file_info() { return file_info_; }

    // Setters
    void set_basic_info(const BasicInfo& info) { basic_info_ = info; }
    void set_space_info(const SpaceInfo& info) { space_info_ = info; }
    void set_mc_info(const MCParams& info) { mc_info_ = info; }
    void set_energy_info(const EnergyInfo& info) { energy_info_ = info; }
    void set_fragment_info(const FragmentInfo& info) { fragment_info_ = info; }
    void set_bias_info(const BiasInfo& info) { bias_info_ = info; }
    void set_file_info(const FileInfo& info) { file_info_ = info; }

    // Convenience methods
    void set_temperature(float temperature) {
        mc_info_.setTemperature(temperature);
    }

    void set_box_size(float x, float y, float z) {
        space_info_.setBoxSize(x, y, z);
    }

    void set_mc_steps(int steps) {
        if (steps <= 0) {
            throw std::invalid_argument("MC steps must be positive");
        }
        mc_info_.mc_steps = steps;
    }

    void set_cutoff(float cutoff) {
        if (cutoff <= 0.0f) {
            throw std::invalid_argument("Cutoff must be positive");
        }
        space_info_.cutoff = cutoff;
        energy_info_.fragment_cutoff = cutoff;
        energy_info_.protein_cutoff = cutoff;
        energy_info_.updateSquaredValues();
    }

    // Update derived values
    void update_derived_values() {
        mc_info_.updateBeta();
        space_info_.updateVolume();
        energy_info_.updateSquaredValues();
        fragment_info_.updateSquaredValues();
        bias_info_.updateSquaredValues();
    }

    // Clear all data
    void clear() {
        basic_info_ = BasicInfo();
        space_info_ = SpaceInfo();
        mc_info_ = MCParams();
        energy_info_ = EnergyInfo();
        fragment_info_ = FragmentInfo();
        bias_info_ = BiasInfo();
        file_info_ = FileInfo();
    }

    // Backward compatibility helper functions
    std::string get_fragment_coordinate_filename(int frag_index) const {
        return file_info_.get_fragment_coordinate_filename(frag_index);
    }

    // String representation
    std::string to_string() const {
        std::stringstream ss;
        ss << "GCMC Parameters:\n";
        ss << "  Version: " << basic_info_.version << "\n";
        ss << "  Temperature: " << mc_info_.temperature << " K\n";
        ss << "  MC Steps: " << mc_info_.mc_steps << "\n";
        ss << "  Box Size: [" << space_info_.box_size[0] << ", " 
           << space_info_.box_size[1] << ", " << space_info_.box_size[2] << "] Å\n";
        ss << "  Cutoff: " << space_info_.cutoff << " Å\n";
        ss << "  Number of fragments: " << file_info_.fragment_names.size();
        return ss.str();
    }

    // Clone method
    std::unique_ptr<Param> clone() const {
        auto cloned = std::make_unique<Param>();
        cloned->basic_info_ = basic_info_;
        cloned->space_info_ = space_info_;
        cloned->mc_info_ = mc_info_;
        cloned->energy_info_ = energy_info_;
        cloned->fragment_info_ = fragment_info_;
        cloned->bias_info_ = bias_info_;
        cloned->file_info_ = file_info_;
        return cloned;
    }

private:
    BasicInfo basic_info_;
    SpaceInfo space_info_;
    MCParams mc_info_;
    EnergyInfo energy_info_;
    FragmentInfo fragment_info_;
    BiasInfo bias_info_;
    FileInfo file_info_;
};

} // namespace param

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_MAIN_HPP 