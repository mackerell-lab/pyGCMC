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
 * @details 包含了GCMC模拟所需的所有参数，包括：
 * - 力场与分子拓扑参数：定义分子间相互作用
 * - 系统几何参数：定义模拟空间
 * - 溶质与溶剂参数：定义化学环境
 * - 算法控制参数：定义采样策略
 * - 模拟控制参数：定义运行条件
 * - 输入输出参数：定义文件路径
 */
class Param {
public:
    // 基本信息结构体
    struct BasicInfo {
        std::string version = "gcmc_v2.0";
        /**
         * @brief 输出详细程度 (verbose)
         * @details 控制程序输出的详细程度：
         * - 0: 只输出基本信息
         * - 1: 输出详细运行信息
         * - 2: 输出调试信息
         */
        int verbosity = 0;                        // verbose, 控制程序输出的详细程度

        /**
         * @brief 调试级别 (debug)
         * @details 控制调试信息的输出级别：
         * - 0: 关闭调试
         * - 1: 基本调试信息
         * - 2: 详细调试信息
         */
        int debug = 0;                           // debug, 控制调试信息的输出级别

        bool print_logfile = false;              // logfile, 是否输出日志文件，用于监控模拟进程
        std::string param_file;                  // paramfile, 参数文件路径，定义模拟的所有参数
        std::string log_file;                    // logfile, 日志文件路径，记录模拟过程
        unsigned int random_seed = 0;            // random_seed, 随机数种子，用于复现性
        int num_threads = 1;                     // nthreads, 并行线程数，用于CPU并行计算
        bool is_box = false;                     // box_size定义标志，定义周期性边界条件范围
        bool init_cycle = false;                 // initcycle, 是否从初始构型开始模拟，如能量最小化后的结构
        bool conserve_fragments = false;         // conserve_frags, 是否固定溶质数量，控制GCMC是否允许动态调整溶质数目
    };

    // 体系空间结构体 - 对应系综参数中的空间参数
    struct SpaceInfo {
        /**
         * @brief 网格间距 (grid_dx)
         * @details 用于空穴偏置采样的网格大小，通常设置为1Å
         * 影响空穴检测的精度和计算效率
         */
        float grid_spacing = 1.0;                // grid_dx, 空穴检测网格分辨率(Å)，用于空穴偏置采样的空间划分

        /**
         * @brief GCMC区域中心坐标
         * @details 定义大正则蒙特卡洛模拟的活性区域中心
         * 用于限制分子插入/删除的空间范围
         */
        std::array<float, 3> gc_center = {0.0, 0.0, 0.0};    // gc_center, GCMC活性区域中心，定义分子插入/删除的空间范围

        std::array<float, 3> sys_center = {0.0, 0.0, 0.0};   // sys_center, 系统中心坐标
        std::array<float, 3> crystal_dim = {0.0, 0.0, 0.0};  // crystal_dim, 晶胞维度参数
        
        /**
         * @brief 模拟盒子大小 (V参数)
         * @details 定义周期性边界条件下的模拟体系大小
         * 对应GCMC系综中的V参数，影响系统的密度和压力
         */
        std::array<float, 3> box_size = {0.0, 0.0, 0.0};     // box_size, 模拟盒子尺寸(Å)，定义周期性边界条件范围

        float volume = 0.0;                      // 计算得到的总体积
        float target_volume = 0.0;               // target_volume, 目标体积
        float sys_box_volume = 0.0;              // 系统盒子体积
        float gcmc_volume = 0.0;                 // GCMC区域体积
        float protein_volume = 0.0;              // 蛋白质体积
        
        /**
         * @brief 网格划分参数
         * @details 控制空穴检测的精确度和效率：
         * - use_vdw_radius_for_grid: 是否使用范德华半径进行网格划分
         * - exclude_hydrogens_from_grid: 是否在网格中排除氢原子
         * - exclude_protein_volume: 是否排除蛋白质体积
         */
        bool use_vdw_radius_for_grid = false;    // use_vdw_radius_for_grid, 是否使用范德华半径进行网格划分
        bool exclude_hydrogens_from_grid = false; // exclude_hydrogens_from_grid, 是否在网格中排除氢原子
        bool exclude_protein_volume = false;      // exclude_protein_volume, 是否排除蛋白质体积
        
        float tmp_prob = 0.0;                    // 临时概率变量
        
        /**
         * @brief 非键相互作用截断距离
         * @details 定义范德华力和静电相互作用的计算截断，通常设为12Å
         * 影响能量计算的精度和效率
         */
        float cutoff = 12.0;                     // cutoff, 非键相互作用截断距离(Å)，影响能量计算精度和效率
    };

    // 蒙特卡洛参数结构体 - 对应运动参数和系综参数
    struct MCInfo {
        /**
         * @brief MC模拟参数
         * @details 控制模拟的基本参数：
         * - mc_steps: 总MC步数
         * - current_step: 当前步数
         * - print_freq: 输出频率
         */
        int mc_steps = 1;                        // mcsteps, MC总步数，若启用系统分区，等效步数=实际步数×分区数
        int current_step = 0;                    // 当前MC步数
        int print_freq = 1;                      // nprint, 输出频率，每隔多少步输出一次日志信息

        /**
         * @brief 热力学参数
         * @details 系统的热力学状态参数：
         * - temperature: 温度T，单位K
         * - beta: β = 1/(kB*T)，用于Metropolis准则
         */
        float temperature = 300.0;               // temperature, 模拟温度(K)，影响Metropolis准则
        float beta = 1.0;                        // 由temperature计算得到，β = 1/(kB*T)

        /**
         * @brief MC移动概率参数
         * @details 控制不同类型MC移动的比例：
         * - insertion_deletion_frac: 插入/删除操作比例
         * - translation_rotation_frac: 平移/旋转操作比例
         */
        float insertion_deletion_frac = 0.5;     // insdel_frac, 插入/删除操作比例
        float translation_rotation_frac = 0.5;    // 平移/旋转操作比例，由insdel_frac计算得到

        /**
         * @brief MC移动参数
         * @details 定义MC移动的最大步长：
         * - max_translation_dist: 最大平移距离，单位Å
         * - max_rotation_angle: 最大旋转角度，单位度
         */
        float max_translation_dist = 1.0;        // 最大平移距离(Å)，控制分子移动步长
        float max_rotation_angle = 30.0;         // 最大旋转角度(度)，控制分子旋转步长

        // MC操作类型列表
        std::vector<std::string> operation_types = {"Ins", "Del", "Trn", "Rot"};  // MC操作类型

        /**
         * @brief MC时间参数
         * @details 控制不同类型移动的时间分配：
         * - mc_time_list: 各类型移动的时间权重
         * - mc_time_cumulative: 累积时间，用于移动类型的选择
         */
        std::vector<float> mc_time_list;         // mctime, MC时间列表，控制不同类型移动的时间分配
        std::vector<float> mc_time_cumulative;   // 累积MC时间，用于移动类型的选择

        /**
         * @brief 操作概率参数
         * @details 不同类型操作的接受概率：
         * - fragment_prob: 片段操作概率，对应Ainsert和Adelete
         * - water_prob: 水分子操作概率
         * - atom_prob: 原子操作概率
         * - test_prob: 测试操作概率
         */
        std::vector<float> fragment_prob;        // 片段操作概率 (对应Ainsert, Adelete)
        std::vector<float> water_prob;           // 水分子操作概率
        std::vector<float> atom_prob;            // 原子操作概率
        std::vector<float> test_prob;            // 测试操作概率

        int rotate_dih_status = 0;               // rotate_dihedral, 二面角旋转状态

        /**
         * @brief 物理常数
         * @details 用于能量计算的常数：
         * - BOLTZMANN: 玻尔兹曼常数kB，单位kcal/mol/K
         * - KCAL_TO_KJ: 能量单位转换因子
         */
        float BOLTZMANN = 0.001987f;             // 玻尔兹曼常数kB (kcal/mol/K)
        float KCAL_TO_KJ = 4.184f;               // 能量单位转换因子 (kcal/mol -> kJ/mol)
    };

    // 能量计算参数结构体 - 对应力场参数
    struct EnergyInfo {
        /**
         * @brief 截断参数
         * @details 控制非键相互作用计算的范围：
         * - use_group_cutoff: 是否使用组截断
         * - fragment_cutoff: 片段能量截断距离
         * - protein_cutoff: 蛋白质能量截断距离
         */
        bool use_group_cutoff = true;            // use_group_cutoff, 是否使用组截断，优化非键相互作用计算
        float fragment_cutoff = 10.0;            // energy_cutoff_frag, 片段能量截断(Å)
        float protein_cutoff = 10.0;             // energy_cutoff_prot, 蛋白质能量截断(Å)
        float fragment_cutoff_squared = 100.0;   // energy_cutoff_frag的平方
        float protein_cutoff_squared = 100.0;    // energy_cutoff_prot的平方

        /**
         * @brief 配对列表参数
         * @details 用于优化非键相互作用计算：
         * - pairlist_cutoff: 配对列表截断距离
         * - pairlist_freq: 配对列表更新频率
         */
        float pairlist_cutoff = 0.0;            // 配对列表截断(Å)，优化非键相互作用计算
        float pairlist_cutoff_squared = 0.0;     // 配对列表截断平方
        unsigned int pairlist_freq = 1000;       // pairlist_freq, 配对列表更新频率

        /**
         * @brief 切换函数参数
         * @details 用于平滑非键相互作用的截断：
         * - use_switching: 是否使用切换函数
         * - switch_dist_fragment: 片段切换距离
         * - switch_dist_protein: 蛋白质切换距离
         */
        bool use_switching = false;              // use_switching, 是否使用切换函数平滑非键相互作用截断
        float switch_dist_fragment = 0.0;        // switch_dist_frag, 片段切换距离(Å)
        float switch_dist_protein = 0.0;         // switch_dist_prot, 蛋白质切换距离(Å)
        float switch_dist_fragment_squared = 0.0; // switch_dist_frag的平方
        float switch_dist_protein_squared = 0.0;  // switch_dist_prot的平方

        /**
         * @brief SW能量参数
         * @details 控制SW能量函数的参数：
         * - energy_sw_ref: SW参考能量
         * - energy_sw_scale: SW能量尺度
         */
        float energy_sw_ref = 1.0;               // SW_reference, SW参考能量
        float energy_sw_scale = 1.0;             // SW_scale, SW能量尺度

        /**
         * @brief SW过滤器参数
         * @details 控制SW过滤器的使用：
         * - test_sw_filters: 是否测试SW过滤器
         * - apply_sw_filters: 是否应用SW过滤器
         */
        bool test_sw_filters = false;            // test_SW_filters, 是否测试SW过滤器
        bool apply_sw_filters = false;           // apply_SW_filters, 是否应用SW过滤器
        bool test_energy = false;                // test_energy, 是否测试能量

        /**
         * @brief 配对列表参数
         * @details 用于优化非键相互作用计算的配对列表参数
         */
        float pair_list_cutoff_fragment = 0.0;   // pair_list_cutoff_frag, 片段配对列表截断(Å)
        float pair_list_cutoff_protein = 0.0;    // pair_list_cutoff_prot, 蛋白质配对列表截断(Å)
        float pair_list_cutoff_fragment_squared = 0.0;  // pair_list_cutoff_frag的平方
        float pair_list_cutoff_protein_squared = 0.0;   // pair_list_cutoff_prot的平方
    };

    // 片段参数结构体 - 对应系综参数和增强振荡μₑₓ协议参数
    struct FragmentInfo {
        /**
         * @brief 水分子参数
         * @details 控制水分子的数量和密度：
         * - water_density: 目标水密度，通常为55.0 M
         * - num_waters: 当前水分子数，对应Ncurrent
         * - target_num_waters: 目标水分子数，对应Ntarget
         */
        float water_density = 55.0;              // 水密度(M)，通常为55.0 M，用于控制水分子数量
        float epsilon = 1.0;                     // epsilon, ε参数
        int num_waters = 0;                      // numwaters, 当前水分子数(对应Ncurrent)
        int target_num_waters = 0;               // target_numwaters, 目标水分子数(对应Ntarget)
        int water_index = 0;                     // sol片段的索引

        /**
         * @brief 过量参数
         * @details 控制分子数量的波动范围：
         * - excess_threshold: 允许的最大偏差比例，对应L参数
         */
        float excess_threshold = 1.0;            // excess_fragments_threshold, 过量阈值(对应L参数)，允许的最大偏差比例

        /**
         * @brief 水分子数均值参数
         * @details 控制水分子数的统计特性：
         * - use_number_water_nbar: 是否使用水分子数均值
         * - use_const_water_nbar: 是否使用固定水分子数均值
         * - const_water_nbar: 固定水分子数均值
         */
        bool use_number_water_nbar = true;       // use_number_water_nbar, 是否使用水分子数均值
        bool use_const_water_nbar = false;       // use_const_water_nbar, 是否使用固定水分子数均值
        int const_water_nbar = 0;                // const_water_nbar, 固定水分子数均值

        /**
         * @brief GCMC截断参数
         * @details 控制GCMC模拟的空间范围：
         * - init_cutoff: 初始化截断距离
         * - gcmc_cutoff: GCMC区域截断距离
         */
        float init_cutoff = 0.0;                 // initial_fragments_cutoff, 初始化截断距离(Å)
        float init_cutoff_squared = 0.0;         // initial_fragments_cutoff的平方
        bool use_gcmc_cutoff = false;            // use_gcmc_cutoff, 是否使用GCMC截断
        float gcmc_cutoff = 0.0;                 // gcmc_cutoff, GCMC区域截断距离(Å)
        float gcmc_cutoff_squared = 0.0;         // gcmc_cutoff的平方

        /**
         * @brief 分子移除列表
         * @details 用于控制分子的移除：
         * - remove_init: 初始化时要移除的分子
         * - remove_excess: 过量时要移除的分子
         */
        std::vector<int> remove_init;            // remove_init, 初始移除列表
        std::vector<int> remove_excess;          // remove_excess, 过量移除列表

        /**
         * @brief 构型和空腔参数
         * @details 用于构型偏置和空穴偏置采样：
         * - confs_list: 可用构型列表
         * - cavity_index_list: 空腔索引列表
         * - cavity_list: 空腔列表
         */
        std::vector<int> confs_list;             // 可用构型列表
        std::vector<int> cavity_index_list;      // 空腔索引列表，用于空穴偏置采样
        std::vector<float> cavity_list;          // 空腔列表，存储空腔大小

        /**
         * @brief 化学势和浓度参数
         * @details 控制GCMC的热力学条件：
         * - conc_list: 目标浓度列表
         * - muex_list: 过量化学势列表，对应μₑₓ
         * - radius_list: 分子半径列表
         */
        std::vector<float> conc_list;            // fragconc, 目标浓度列表(M)，如溶质(0.25 M)和溶剂(55 M)
        std::vector<float> muex_list;            // fragmuex, 过量化学势列表(kcal/mol)，用于控制插入/删除概率
        std::vector<float> radius_list;          // fragradius, 分子半径列表(Å)

        std::vector<int> conf_list;              // 构型列表，用于构型偏置采样
        int flag_remove_init = 0;                // remove_init标志，初始移除标志
        int flag_remove_excess = 0;              // remove_excess标志，过量移除标志
        int total_protitp_size = 0;              // protitp文件总数，蛋白质拓扑文件总数
        std::vector<int> fragconf_list;          // fragconfs, 片段构型列表
    };

    // 偏置采样参数结构体 - 对应偏置采样参数
    struct BiasInfo {
        /**
         * @brief 空穴偏置参数
         * @details 用于提高插入/删除移动的接受率：
         * - use_cavity_bias: 是否使用空穴偏置采样
         * - sigma: 空穴大小参数
         */
        bool use_cavity_bias = false;            // use_cavity_bias, 是否使用空穴偏置采样，仅向空腔区域尝试插入以提升接受率
        float sigma = 2.4;                       // sigma, σ参数，空穴大小参数(Å)
        float sigma_squared = 5.76;              // sigma的平方，σ²参数

        /**
         * @brief 构型偏置参数
         * @details 用于提高构型采样效率：
         * - use_conf_bias: 是否使用构型偏置采样
         * - num_conf_bias_trials: 每次尝试的构型数，对应n参数
         */
        bool use_conf_bias = false;              // use_conf_bias, 是否使用构型偏置采样，每次插入尝试多构型并按能量权重选择
        unsigned int num_conf_bias_trials = 10;   // num_conf_bias_trial, 构型偏置尝试次数(对应n参数)
    };

    // 文件路径参数结构体 - 对应模拟控制参数
    struct FileInfo {
        /**
         * @brief 主要输入文件
         * @details 模拟所需的基本文件：
         * - topology_file: 系统拓扑文件
         * - input_pdb_file: 初始构型文件
         * - output_pdb_file: 轨迹输出文件
         */
        std::string topology_file;               // top, 系统拓扑文件，包含分子类型、原子列表、力场参数
        std::string input_pdb_file;              // pdb, 初始构型文件，定义原子坐标和分子排布
        std::string output_pdb_file;             // op_pdb, 轨迹输出文件，保存模拟后的构型
        std::string output_top_file;             // op_top, 输出拓扑文件

        /**
         * @brief 力场文件
         * @details 力场参数和原子类型文件：
         * - atomtype_file: 原子类型定义文件
         * - par_files: 力场参数文件列表
         */
        std::string atomtype_file;               // atomtypes, 原子类型定义文件，映射原子名称到力场类型
        std::string monomer_dir;                 // monomerdir, 单体目录
        std::string conc_norm = "water";         // conc_norm, 浓度归一化方式
        std::string conc_region = "total";       // conc_region, 浓度计算区域
        std::vector<std::string> par_files;      // par, 力场参数文件，定义非键相互作用参数

        /**
         * @brief 拓扑文件
         * @details 不同组分的拓扑文件：
         * - protein_top_files: 蛋白质拓扑文件
         * - fragment_top_files: 片段拓扑文件
         * - fragment_names: 片段名称列表
         */
        std::vector<std::string> protein_top_files;  // protitp, 蛋白质拓扑文件
        std::vector<std::string> fragment_top_files; // fragitp, 溶质分子拓扑文件，包含原子类型、电荷、键合参数
        std::vector<std::string> fragment_names;     // fragname, 溶质和溶剂名称列表
        std::vector<std::string> fragment_mqtr_files;// fragmqtr, 片段MQTR文件
        std::string tmp_frag_name;               // 临时片段名称

        /**
         * @brief 映射文件参数
         * @details 用于生成和存储空间映射：
         * - generate_maps: 是否生成映射文件
         * - map_prefix: 映射文件名前缀
         */
        bool generate_maps = false;              // map_generation, 是否生成空间映射文件
        std::string map_prefix = "gc_maps";      // map_filename_prefix, 映射文件名前缀
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

 