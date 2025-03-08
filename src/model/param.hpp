// src/model/param.hpp

#pragma once
#ifndef PYGCMC_MODEL_PARAM_HPP
#define PYGCMC_MODEL_PARAM_HPP

#include <vector>
#include <string>
#include <array>
#include <memory>
#include <stdexcept>  // for std::runtime_error
#include <cmath>      // for std::sqrt
#include <algorithm>  // for std::find, std::lower_bound, std::distance

namespace pygcmc {
namespace model {

/**
 * @brief GCMC parameter class for storing all parameters needed for calculation
 * @details Contains all parameters required for GCMC simulation, including:
 * - Force field and molecular topology parameters: define intermolecular interactions
 * - System geometry parameters: define simulation space
 * - Solute and solvent parameters: define chemical environment
 * - Algorithm control parameters: define sampling strategy
 * - Simulation control parameters: define running conditions
 * - Input/output parameters: define file paths
 */
class Param {
public:
    // Basic information structure
    struct BasicInfo {
        std::string version = "gcmc_v2.0";
        /**
         * @brief Output verbosity level
         * @details Controls the level of detail in program output:
         * - 0: Output only basic information
         * - 1: Output detailed running information
         * - 2: Output debugging information
         */
        int verbosity = 0;                        // verbose, controls the level of detail in program output

        /**
         * @brief Debug level
         * @details Controls the output level of debugging information:
         * - 0: Disable debugging
         * - 1: Basic debugging information
         * - 2: Detailed debugging information
         */
        int debug = 0;                           // debug, controls the output level of debugging information

        bool print_logfile = false;              // logfile, whether to output log file for monitoring simulation process
        std::string param_file;                  // paramfile, parameter file path, defines all simulation parameters
        std::string log_file;                    // logfile, log file path, records simulation process
        unsigned int random_seed = 0;            // random_seed, random number seed for reproducibility
        int num_threads = 1;                     // nthreads, number of parallel threads for CPU parallel computation
        bool is_box = false;                     // box_size flag, defines periodic boundary condition range
        bool init_cycle = false;                 // initcycle, whether to start simulation from initial configuration, e.g., structure after energy minimization
        bool conserve_fragments = false;         // conserve_frags, whether to fix solute quantity, controls if GCMC allows dynamic adjustment of solute number
    };

    // System space structure - corresponds to spatial parameters in ensemble parameters
    struct SpaceInfo {
        /**
         * @brief Grid spacing (grid_dx)
         * @details Grid size for cavity-biased sampling, typically set to 1Å
         * Affects the precision and computational efficiency of cavity detection
         */
        float grid_spacing = 1.0;                // grid_dx, cavity detection grid resolution (Å), used for spatial division in cavity-biased sampling

        /**
         * @brief GCMC region center coordinates
         * @details Defines the center of the active region for grand canonical Monte Carlo simulation
         * Used to limit the spatial range for molecule insertion/deletion
         */
        std::array<float, 3> gc_center = {0.0, 0.0, 0.0};    // gc_center, GCMC active region center, defines spatial range for molecule insertion/deletion

        std::array<float, 3> sys_center = {0.0, 0.0, 0.0};   // sys_center, system center coordinates
        std::array<float, 3> crystal_dim = {0.0, 0.0, 0.0};  // crystal_dim, crystal cell dimension parameters
        
        /**
         * @brief Simulation box size (V parameter)
         * @details Defines the simulation system size under periodic boundary conditions
         * Corresponds to the V parameter in GCMC ensemble, affects system density and pressure
         */
        std::array<float, 3> box_size = {0.0, 0.0, 0.0};     // box_size, simulation box dimensions (Å), defines periodic boundary condition range

        float volume = 0.0;                      // calculated total volume
        float target_volume = 0.0;               // target_volume, target volume
        float sys_box_volume = 0.0;              // system box volume
        float gcmc_volume = 0.0;                 // GCMC region volume
        float protein_volume = 0.0;              // protein volume
        
        /**
         * @brief Grid division parameters
         * @details Control the precision and efficiency of cavity detection:
         * - use_vdw_radius_for_grid: whether to use van der Waals radius for grid division
         * - exclude_hydrogens_from_grid: whether to exclude hydrogen atoms from the grid
         * - exclude_protein_volume: whether to exclude protein volume
         */
        bool use_vdw_radius_for_grid = false;    // use_vdw_radius_for_grid, whether to use van der Waals radius for grid division
        bool exclude_hydrogens_from_grid = false; // exclude_hydrogens_from_grid, whether to exclude hydrogen atoms from the grid
        bool exclude_protein_volume = false;      // exclude_protein_volume, whether to exclude protein volume
        
        float tmp_prob = 0.0;                    // temporary probability variable
        
        /**
         * @brief Non-bonded interaction cutoff distance
         * @details Defines the calculation cutoff for van der Waals and electrostatic interactions, typically set to 12Å
         * Affects the precision and efficiency of energy calculation
         */
        float cutoff = 12.0;                     // cutoff, non-bonded interaction cutoff distance (Å), affects energy calculation precision and efficiency
    };

    // Monte Carlo parameter structure - corresponds to movement parameters and ensemble parameters
    struct MCInfo {
        /**
         * @brief MC simulation parameters
         * @details Control basic simulation parameters:
         * - mc_steps: total MC steps
         * - current_step: current step
         * - print_freq: output frequency
         */
        int mc_steps = 1;                        // mcsteps, total MC steps, if system partitioning is enabled, equivalent steps = actual steps × number of partitions
        int current_step = 0;                    // current MC step
        int print_freq = 1;                      // nprint, output frequency, output log information every n steps

        /**
         * @brief Thermodynamic parameters
         * @details Thermodynamic state parameters of the system:
         * - temperature: temperature T, unit K
         * - beta: β = 1/(kB*T), used for Metropolis criterion
         */
        float temperature = 300.0;               // temperature, simulation temperature (K), affects Metropolis criterion
        float beta = 1.0;                        // calculated from temperature, β = 1/(kB*T)

        /**
         * @brief MC move probability parameters
         * @details Control the ratio of different types of MC moves:
         * - insertion_deletion_frac: insertion/deletion operation ratio
         * - translation_rotation_frac: translation/rotation operation ratio
         */
        float insertion_deletion_frac = 0.5;     // insdel_frac, insertion/deletion operation ratio
        float translation_rotation_frac = 0.5;    // translation/rotation operation ratio, calculated from insdel_frac

        /**
         * @brief MC move parameters
         * @details Define the maximum step size for MC moves:
         * - max_translation_dist: maximum translation distance, unit Å
         * - max_rotation_angle: maximum rotation angle, unit degrees
         */
        float max_translation_dist = 1.0;        // maximum translation distance (Å), controls molecule movement step size
        float max_rotation_angle = 30.0;         // maximum rotation angle (degrees), controls molecule rotation step size

        // MC operation type list
        std::vector<std::string> operation_types = {"Ins", "Del", "Trn", "Rot"};  // MC operation types

        /**
         * @brief MC time parameters
         * @details Control time allocation for different types of moves:
         * - mc_time_list: time weights for each type of move
         * - mc_time_cumulative: cumulative time, used for move type selection
         */
        std::vector<float> mc_time_list;         // mctime, MC time list, controls time allocation for different types of moves
        std::vector<float> mc_time_cumulative;   // cumulative MC time, used for move type selection

        /**
         * @brief Operation probability parameters
         * @details Acceptance probabilities for different types of operations:
         * - fragment_prob: fragment operation probability, corresponding to Ainsert and Adelete
         * - water_prob: water molecule operation probability
         * - atom_prob: atom operation probability
         * - test_prob: test operation probability
         */
        std::vector<float> fragment_prob;        // fragment operation probability (corresponding to Ainsert, Adelete)
        std::vector<float> water_prob;           // water molecule operation probability
        std::vector<float> atom_prob;            // atom operation probability
        std::vector<float> test_prob;            // test operation probability

        int rotate_dih_status = 0;               // rotate_dihedral, dihedral rotation status

        /**
         * @brief Physical constants
         * @details Constants used for energy calculation:
         * - BOLTZMANN: Boltzmann constant kB, unit kcal/mol/K
         * - KCAL_TO_KJ: energy unit conversion factor
         */
        float BOLTZMANN = 0.001987f;             // Boltzmann constant kB (kcal/mol/K)
        float KCAL_TO_KJ = 4.184f;               // energy unit conversion factor (kcal/mol -> kJ/mol)
    };

    // Energy calculation parameter structure - corresponds to force field parameters
    struct EnergyInfo {
        /**
         * @brief Cutoff parameters
         * @details Control the range of non-bonded interaction calculations:
         * - use_group_cutoff: whether to use group cutoff
         * - fragment_cutoff: fragment energy cutoff distance
         * - protein_cutoff: protein energy cutoff distance
         */
        bool use_group_cutoff = true;            // use_group_cutoff, whether to use group cutoff, optimize non-bonded interaction calculation
        float fragment_cutoff = 10.0;            // energy_cutoff_frag, fragment energy cutoff (Å)
        float protein_cutoff = 10.0;             // energy_cutoff_prot, protein energy cutoff (Å)
        float fragment_cutoff_squared = 100.0;   // square of energy_cutoff_frag
        float protein_cutoff_squared = 100.0;    // square of energy_cutoff_prot

        /**
         * @brief Pairlist parameters
         * @details Used to optimize non-bonded interaction calculations:
         * - pairlist_cutoff: pairlist cutoff distance
         * - pairlist_freq: pairlist update frequency
         */
        float pairlist_cutoff = 0.0;            // pairlist cutoff (Å), optimize non-bonded interaction calculation
        float pairlist_cutoff_squared = 0.0;     // square of pairlist cutoff
        unsigned int pairlist_freq = 1000;       // pairlist_freq, pairlist update frequency

        /**
         * @brief Switching function parameters
         * @details Used to smooth non-bonded interaction cutoff:
         * - use_switching: whether to use switching function
         * - switch_dist_fragment: fragment switching distance
         * - switch_dist_protein: protein switching distance
         */
        bool use_switching = false;              // use_switching, whether to use switching function to smooth non-bonded interaction cutoff
        float switch_dist_fragment = 0.0;        // switch_dist_frag, fragment switching distance (Å)
        float switch_dist_protein = 0.0;         // switch_dist_prot, protein switching distance (Å)
        float switch_dist_fragment_squared = 0.0; // square of switch_dist_frag
        float switch_dist_protein_squared = 0.0;  // square of switch_dist_prot

        /**
         * @brief SW energy parameters
         * @details Control SW energy function parameters:
         * - energy_sw_ref: SW reference energy
         * - energy_sw_scale: SW energy scale
         */
        float energy_sw_ref = 1.0;               // SW_reference, SW reference energy
        float energy_sw_scale = 1.0;             // SW_scale, SW energy scale

        /**
         * @brief SW filter parameters
         * @details Control the use of SW filter:
         * - test_sw_filters: whether to test SW filters
         * - apply_sw_filters: whether to apply SW filters
         */
        bool test_sw_filters = false;            // test_SW_filters, whether to test SW filters
        bool apply_sw_filters = false;           // apply_SW_filters, whether to apply SW filters
        bool test_energy = false;                // test_energy, whether to test energy

        /**
         * @brief Pairlist parameters
         * @details Pairlist parameters for optimizing non-bonded interaction calculations
         */
        float pair_list_cutoff_fragment = 0.0;   // pair_list_cutoff_frag, fragment pairlist cutoff (Å)
        float pair_list_cutoff_protein = 0.0;    // pair_list_cutoff_prot, protein pairlist cutoff (Å)
        float pair_list_cutoff_fragment_squared = 0.0;  // square of pair_list_cutoff_frag
        float pair_list_cutoff_protein_squared = 0.0;   // square of pair_list_cutoff_prot
    };

    // Fragment parameter structure - corresponds to ensemble parameters and enhanced oscillation μₑₓ protocol parameters
    struct FragmentInfo {
        /**
         * @brief Water molecule parameters
         * @details Control water molecule quantity and density:
         * - water_density: target water density, typically 55.0 M
         * - num_waters: current number of water molecules, corresponding to Ncurrent
         * - target_num_waters: target number of water molecules, corresponding to Ntarget
         */
        float water_density = 55.0;              // water density (M), typically 55.0 M, used to control water molecule quantity
        float epsilon = 1.0;                     // epsilon, ε parameter
        int num_waters = 0;                      // numwaters, current number of water molecules (corresponding to Ncurrent)
        int target_num_waters = 0;               // target_numwaters, target number of water molecules (corresponding to Ntarget)
        int water_index = 0;                     // index of sol fragment

        /**
         * @brief Excess parameters
         * @details Control the fluctuation range of molecule quantity:
         * - excess_threshold: maximum allowed deviation ratio, corresponding to L parameter
         */
        float excess_threshold = 1.0;            // excess_fragments_threshold, excess threshold (corresponding to L parameter), maximum allowed deviation ratio

        /**
         * @brief Water molecule number mean parameters
         * @details Control the statistical characteristics of water molecule number:
         * - use_number_water_nbar: whether to use water molecule number mean
         * - use_const_water_nbar: whether to use constant water molecule number mean
         * - const_water_nbar: constant water molecule number mean
         */
        bool use_number_water_nbar = true;       // use_number_water_nbar, whether to use water molecule number mean
        bool use_const_water_nbar = false;       // use_const_water_nbar, whether to use constant water molecule number mean
        int const_water_nbar = 0;                // const_water_nbar, constant water molecule number mean

        /**
         * @brief GCMC cutoff parameters
         * @details Control the spatial range of GCMC simulation:
         * - init_cutoff: initialization cutoff distance
         * - gcmc_cutoff: GCMC region cutoff distance
         */
        float init_cutoff = 0.0;                 // initial_fragments_cutoff, initialization cutoff distance (Å)
        float init_cutoff_squared = 0.0;         // square of initial_fragments_cutoff
        bool use_gcmc_cutoff = false;            // use_gcmc_cutoff, whether to use GCMC cutoff
        float gcmc_cutoff = 0.0;                 // gcmc_cutoff, GCMC region cutoff distance (Å)
        float gcmc_cutoff_squared = 0.0;         // square of gcmc_cutoff

        /**
         * @brief Molecule removal list
         * @details Used to control molecule removal:
         * - remove_init: molecules to remove during initialization
         * - remove_excess: molecules to remove when in excess
         */
        std::vector<int> remove_init;            // remove_init, initial removal list
        std::vector<int> remove_excess;          // remove_excess, excess removal list

        /**
         * @brief Configuration and cavity parameters
         * @details Used for configuration bias and cavity-biased sampling:
         * - confs_list: available configuration list
         * - cavity_index_list: cavity index list
         * - cavity_list: cavity list
         */
        std::vector<int> confs_list;             // available configuration list
        std::vector<int> cavity_index_list;      // cavity index list, used for cavity-biased sampling
        std::vector<float> cavity_list;          // cavity list, stores cavity sizes

        /**
         * @brief Chemical potential and concentration parameters
         * @details Control thermodynamic conditions for GCMC:
         * - conc_list: target concentration list
         * - muex_list: excess chemical potential list, corresponding to μₑₓ
         * - radius_list: molecule radius list
         */
        std::vector<float> conc_list;            // fragconc, target concentration list (M), such as solute (0.25 M) and solvent (55 M)
        std::vector<float> muex_list;            // fragmuex, excess chemical potential list (kcal/mol), used to control insertion/deletion probability
        std::vector<float> radius_list;          // fragradius, molecule radius list (Å)

        std::vector<int> conf_list;              // configuration list, used for configuration bias sampling
        int flag_remove_init = 0;                // remove_init flag, initial removal flag
        int flag_remove_excess = 0;              // remove_excess flag, excess removal flag
        int total_protitp_size = 0;              // total number of protitp files, total number of protein topology files
        std::vector<int> fragconf_list;          // fragconfs, fragment configuration list
    };

    // Biased sampling parameter structure - corresponds to biased sampling parameters
    struct BiasInfo {
        /**
         * @brief Cavity bias parameters
         * @details Used to improve the acceptance rate of insertion/deletion moves:
         * - use_cavity_bias: whether to use cavity-biased sampling
         * - sigma: cavity size parameter
         */
        bool use_cavity_bias = false;            // use_cavity_bias, whether to use cavity-biased sampling, only attempt insertions into cavity regions to improve acceptance rate
        float sigma = 2.4;                       // sigma, σ parameter, cavity size parameter (Å)
        float sigma_squared = 5.76;              // square of sigma, σ² parameter

        /**
         * @brief Configuration bias parameters
         * @details Used to improve configuration sampling efficiency:
         * - use_conf_bias: whether to use configuration-biased sampling
         * - num_conf_bias_trials: number of configurations to try each time, corresponding to n parameter
         */
        bool use_conf_bias = false;              // use_conf_bias, whether to use configuration-biased sampling, try multiple configurations for each insertion and select based on energy weights
        unsigned int num_conf_bias_trials = 10;   // num_conf_bias_trial, configuration bias trial count (corresponding to n parameter)
    };

    // File path parameter structure - corresponds to simulation control parameters
    struct FileInfo {
        /**
         * @brief Main input files
         * @details Basic files required for simulation:
         * - topology_file: system topology file
         * - input_pdb_file: initial configuration file
         * - output_pdb_file: trajectory output file
         */
        std::string topology_file;               // top, system topology file, containing molecule types, atom lists, force field parameters
        std::string input_pdb_file;              // pdb, initial configuration file, defining atom coordinates and molecule arrangement
        std::string output_pdb_file;             // op_pdb, trajectory output file, saving post-simulation configurations
        std::string output_top_file;             // op_top, output topology file

        /**
         * @brief Force field files
         * @details Force field parameters and atom type files:
         * - atomtype_file: atom type definition file
         * - par_files: force field parameter file list
         */
        std::string atomtype_file;               // atomtypes, atom type definition file, mapping atom names to force field types
        std::string monomer_dir;                 // monomerdir, monomer directory
        std::string conc_norm = "water";         // conc_norm, concentration normalization method
        std::string conc_region = "total";       // conc_region, concentration calculation region
        std::vector<std::string> par_files;      // par, force field parameter files, defining non-bonded interaction parameters

        /**
         * @brief Topology files
         * @details Topology files for different components:
         * - protein_top_files: protein topology files
         * - fragment_top_files: fragment topology files
         * - fragment_names: fragment name list
         */
        std::vector<std::string> protein_top_files;  // protitp, protein topology files
        std::vector<std::string> fragment_top_files; // fragitp, solute molecule topology files, containing atom types, charges, bonded parameters
        std::vector<std::string> fragment_names;     // fragname, solute and solvent name list
        std::vector<std::string> fragment_mqtr_files;// fragmqtr, fragment MQTR files
        std::string tmp_frag_name;               // temporary fragment name

        /**
         * @brief Mapping file parameters
         * @details Used to generate and store spatial mapping:
         * - generate_maps: whether to generate mapping files
         * - map_prefix: mapping file name prefix
         */
        bool generate_maps = false;              // map_generation, whether to generate spatial mapping files
        std::string map_prefix = "gc_maps";      // map_filename_prefix, mapping file name prefix
    };

    // Constructor
    Param() = default;

    // Getters
    const BasicInfo& get_basic_info() const { return basic_info_; }
    const SpaceInfo& get_space_info() const { return space_info_; }
    const MCInfo& get_mc_info() const { return mc_info_; }
    const EnergyInfo& get_energy_info() const { return energy_info_; }
    const FragmentInfo& get_fragment_info() const { return fragment_info_; }
    const BiasInfo& get_bias_info() const { return bias_info_; }
    const FileInfo& get_file_info() const { return file_info_; }

    // Non-const getters for modification
    BasicInfo& get_basic_info() { return basic_info_; }
    SpaceInfo& get_space_info() { return space_info_; }
    MCInfo& get_mc_info() { return mc_info_; }
    EnergyInfo& get_energy_info() { return energy_info_; }
    FragmentInfo& get_fragment_info() { return fragment_info_; }
    BiasInfo& get_bias_info() { return bias_info_; }
    FileInfo& get_file_info() { return file_info_; }

    // Setters
    void set_basic_info(const BasicInfo& info) { basic_info_ = info; }
    void set_space_info(const SpaceInfo& info) { space_info_ = info; }
    void set_mc_info(const MCInfo& info) { mc_info_ = info; }
    void set_energy_info(const EnergyInfo& info) { energy_info_ = info; }
    void set_fragment_info(const FragmentInfo& info) { fragment_info_ = info; }
    void set_bias_info(const BiasInfo& info) { bias_info_ = info; }
    void set_file_info(const FileInfo& info) { file_info_ = info; }

    // Clear all data
    void clear() {
        basic_info_ = BasicInfo();
        space_info_ = SpaceInfo();
        mc_info_ = MCInfo();
        energy_info_ = EnergyInfo();
        fragment_info_ = FragmentInfo();
        bias_info_ = BiasInfo();
        file_info_ = FileInfo();
    }

    // Simple helper functions
    std::string get_fragment_coordinate_filename(int frag_index) const {
        return file_info_.monomer_dir + "/" + file_info_.fragment_names[frag_index] + ".pdb";
    }

private:
    BasicInfo basic_info_;
    SpaceInfo space_info_;
    MCInfo mc_info_;
    EnergyInfo energy_info_;
    FragmentInfo fragment_info_;
    BiasInfo bias_info_;
    FileInfo file_info_;
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_HPP

 