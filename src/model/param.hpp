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
 * @brief GCMC参数类，用于存储所有计算所需的参数
 */
class Param {
public:
    // 基本信息结构体
    struct BasicInfo {
        std::string version = "gcmc_v2.0";
        int verbosity = 0;                        // verbose
        int debug = 0;                           // debug
        bool print_logfile = false;              // logfile
        std::string param_file;                  // paramfile
        std::string log_file;                    // logfile
        unsigned int random_seed = 0;            // random_seed
        int num_threads = 1;                     // nthreads
        bool is_box = false;                     // box_size定义标志
        bool init_cycle = false;                 // initcycle
        bool conserve_fragments = false;         // conserve_frags
    };

    // 体系空间结构体
    struct SpaceInfo {
        float grid_spacing = 1.0;                // grid_dx
        std::array<float, 3> gc_center = {0.0, 0.0, 0.0};    // gc_center
        std::array<float, 3> sys_center = {0.0, 0.0, 0.0};   // sys_center
        std::array<float, 3> crystal_dim = {0.0, 0.0, 0.0};  // crystal_dim
        std::array<float, 3> box_size = {0.0, 0.0, 0.0};     // box_size
        float volume = 0.0;                      // 计算得到的体积
        float target_volume = 0.0;               // target_volume
        float sys_box_volume = 0.0;              // 系统盒子体积
        float gcmc_volume = 0.0;                 // GCMC区域体积
        float protein_volume = 0.0;              // 蛋白质体积
        bool use_vdw_radius_for_grid = false;    // use_vdw_radius_for_grid
        bool exclude_hydrogens_from_grid = false; // exclude_hydrogens_from_grid
        bool exclude_protein_volume = false;      // exclude_protein_volume
        float tmp_prob = 0.0;                    // 临时概率变量
        float cutoff = 12.0;                     // cutoff
    };

    // 蒙特卡洛参数结构体
    struct MCInfo {
        int mc_steps = 1;                        // mcsteps
        int current_step = 0;                    // 当前步数
        int print_freq = 1;                      // nprint
        float temperature = 300.0;               // temperature
        float beta = 1.0;                        // 由temperature计算得到
        float insertion_deletion_frac = 0.5;     // insdel_frac
        float translation_rotation_frac = 0.5;    // 由insdel_frac计算得到
        float max_translation_dist = 1.0;        // 最大平移距离
        float max_rotation_angle = 30.0;         // 最大旋转角度
        std::vector<std::string> operation_types = {"Ins", "Del", "Trn", "Rot"};  // MC操作类型
        std::vector<float> mc_time_list;         // mctime
        std::vector<float> mc_time_cumulative;   // 累积mctime
        std::vector<float> fragment_prob;        // 片段操作概率
        std::vector<float> water_prob;           // 水分子操作概率
        std::vector<float> atom_prob;            // 原子操作概率
        std::vector<float> test_prob;            // 测试操作概率
        int rotate_dih_status = 0;               // rotate_dihedral
        float BOLTZMANN = 0.001987f;             // 玻尔兹曼常数
        float KCAL_TO_KJ = 4.184f;               // 能量单位转换因子
    };

    // 能量计算参数结构体
    struct EnergyInfo {
        bool use_group_cutoff = true;            // use_group_cutoff
        float fragment_cutoff = 10.0;            // energy_cutoff_frag
        float protein_cutoff = 10.0;             // energy_cutoff_prot
        float fragment_cutoff_squared = 100.0;   // energy_cutoff_frag的平方
        float protein_cutoff_squared = 100.0;    // energy_cutoff_prot的平方
        float pairlist_cutoff = 0.0;            // 配对列表截断
        float pairlist_cutoff_squared = 0.0;     // 配对列表截断平方
        unsigned int pairlist_freq = 1000;       // pairlist_freq
        bool use_switching = false;              // use_switching
        float switch_dist_fragment = 0.0;        // switch_dist_frag
        float switch_dist_protein = 0.0;         // switch_dist_prot
        float switch_dist_fragment_squared = 0.0; // switch_dist_frag的平方
        float switch_dist_protein_squared = 0.0;  // switch_dist_prot的平方
        float energy_sw_ref = 1.0;               // SW_reference
        float energy_sw_scale = 1.0;             // SW_scale
        bool test_sw_filters = false;            // test_SW_filters
        bool apply_sw_filters = false;           // apply_SW_filters
        bool test_energy = false;                // test_energy
        float pair_list_cutoff_fragment = 0.0;   // pair_list_cutoff_frag
        float pair_list_cutoff_protein = 0.0;    // pair_list_cutoff_prot
        float pair_list_cutoff_fragment_squared = 0.0;  // pair_list_cutoff_frag的平方
        float pair_list_cutoff_protein_squared = 0.0;   // pair_list_cutoff_prot的平方
    };

    // 片段参数结构体
    struct FragmentInfo {
        float water_density = 55.0;              // sol片段的fragconc值
        float epsilon = 1.0;                     // epsilon
        int num_waters = 0;                      // numwaters
        int target_num_waters = 0;               // target_numwaters
        int water_index = 0;                     // sol片段的索引
        float excess_threshold = 1.0;            // excess_fragments_threshold
        bool use_number_water_nbar = true;       // use_number_water_nbar
        bool use_const_water_nbar = false;       // use_const_water_nbar
        int const_water_nbar = 0;                // const_water_nbar
        float init_cutoff = 0.0;                 // initial_fragments_cutoff
        float init_cutoff_squared = 0.0;         // initial_fragments_cutoff的平方
        bool use_gcmc_cutoff = false;            // use_gcmc_cutoff
        float gcmc_cutoff = 0.0;                 // gcmc_cutoff
        float gcmc_cutoff_squared = 0.0;         // gcmc_cutoff的平方
        std::vector<int> remove_init;            // remove_init
        std::vector<int> remove_excess;          // remove_excess
        std::vector<int> confs_list;             // 构型列表
        std::vector<int> cavity_index_list;      // 空腔索引列表
        std::vector<float> cavity_list;          // 空腔列表
        std::vector<float> conc_list;            // fragconc
        std::vector<float> muex_list;            // fragmuex
        std::vector<float> radius_list;          // fragradius
        std::vector<int> conf_list;              // 构型列表
        int flag_remove_init = 0;                // remove_init标志
        int flag_remove_excess = 0;              // remove_excess标志
        int total_protitp_size = 0;              // protitp文件总数
        std::vector<int> fragconf_list;          // fragconfs
    };

    // 偏置采样参数结构体
    struct BiasInfo {
        bool use_cavity_bias = false;            // use_cavity_bias
        float sigma = 2.4;                       // sigma
        float sigma_squared = 5.76;              // sigma的平方
        bool use_conf_bias = false;              // use_conf_bias
        unsigned int num_conf_bias_trials = 10;   // num_conf_bias_trial
    };

    // 文件路径参数结构体
    struct FileInfo {
        std::string topology_file;               // top
        std::string input_pdb_file;              // pdb
        std::string output_pdb_file;             // op_pdb
        std::string output_top_file;             // op_top
        std::string atomtype_file;               // atomtypes
        std::string monomer_dir;                 // monomerdir
        std::string conc_norm = "water";         // conc_norm
        std::string conc_region = "total";       // conc_region
        std::vector<std::string> par_files;      // par
        std::vector<std::string> protein_top_files;  // protitp
        std::vector<std::string> fragment_top_files; // fragitp
        std::vector<std::string> fragment_names;     // fragname
        std::vector<std::string> fragment_mqtr_files;// fragmqtr
        std::string tmp_frag_name;               // 临时片段名称
        bool generate_maps = false;              // map_generation
        std::string map_prefix = "gc_maps";      // map_filename_prefix
    };

    // 构造函数
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

    // 清除所有数据
    void clear() {
        basic_info_ = BasicInfo();
        space_info_ = SpaceInfo();
        mc_info_ = MCInfo();
        energy_info_ = EnergyInfo();
        fragment_info_ = FragmentInfo();
        bias_info_ = BiasInfo();
        file_info_ = FileInfo();
    }

    // 简单的辅助函数
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

 