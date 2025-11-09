// src/bindings/model/ModelParamBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;
using namespace pygcmc::model;

namespace pygcmc {
namespace bindings {
namespace model {

void init_param_bindings(py::module&, py::module& model_module) {
    // Bind Param class and its nested structs
    auto param = py::class_<::pygcmc::model::Param>(model_module, "Param")
        .def(py::init<>())
        // Basic info
        .def_property("basic_info",
            py::overload_cast<>(&::pygcmc::model::Param::get_basic_info, py::const_),
            py::overload_cast<>(&::pygcmc::model::Param::get_basic_info))
        // Space info
        .def_property("space_info",
            py::overload_cast<>(&::pygcmc::model::Param::get_space_info, py::const_),
            py::overload_cast<>(&::pygcmc::model::Param::get_space_info))
        // MC info
        .def_property("mc_info",
            py::overload_cast<>(&::pygcmc::model::Param::get_mc_info, py::const_),
            py::overload_cast<>(&::pygcmc::model::Param::get_mc_info))
        // Energy info
        .def_property("energy_info",
            py::overload_cast<>(&::pygcmc::model::Param::get_energy_info, py::const_),
            py::overload_cast<>(&::pygcmc::model::Param::get_energy_info))
        // Fragment info
        .def_property("fragment_info",
            py::overload_cast<>(&::pygcmc::model::Param::get_fragment_info, py::const_),
            py::overload_cast<>(&::pygcmc::model::Param::get_fragment_info))
        // Bias info
        .def_property("bias_info",
            py::overload_cast<>(&::pygcmc::model::Param::get_bias_info, py::const_),
            py::overload_cast<>(&::pygcmc::model::Param::get_bias_info))
        // File info
        .def_property("file_info",
            py::overload_cast<>(&::pygcmc::model::Param::get_file_info, py::const_),
            py::overload_cast<>(&::pygcmc::model::Param::get_file_info))
        .def("clear", &::pygcmc::model::Param::clear);

    // Bind BasicInfo struct
    py::class_<::pygcmc::model::Param::BasicInfo>(param, "BasicInfo")
        .def(py::init<>())
        .def_readwrite("version", &::pygcmc::model::Param::BasicInfo::version)
        .def_readwrite("verbosity", &::pygcmc::model::Param::BasicInfo::verbosity)
        .def_readwrite("debug", &::pygcmc::model::Param::BasicInfo::debug)
        .def_readwrite("print_logfile", &::pygcmc::model::Param::BasicInfo::print_logfile)
        .def_readwrite("param_file", &::pygcmc::model::Param::BasicInfo::param_file)
        .def_readwrite("log_file", &::pygcmc::model::Param::BasicInfo::log_file)
        .def_readwrite("random_seed", &::pygcmc::model::Param::BasicInfo::random_seed)
        .def_readwrite("num_threads", &::pygcmc::model::Param::BasicInfo::num_threads)
        .def_readwrite("is_box", &::pygcmc::model::Param::BasicInfo::is_box)
        .def_readwrite("init_cycle", &::pygcmc::model::Param::BasicInfo::init_cycle)
        .def_readwrite("conserve_fragments", &::pygcmc::model::Param::BasicInfo::conserve_fragments);

    // Bind SpaceInfo struct
    py::class_<::pygcmc::model::Param::SpaceInfo>(param, "SpaceInfo")
        .def(py::init<>())
        .def_readwrite("grid_spacing", &::pygcmc::model::Param::SpaceInfo::grid_spacing)
        .def_readwrite("gc_center", &::pygcmc::model::Param::SpaceInfo::gc_center)
        .def_readwrite("sys_center", &::pygcmc::model::Param::SpaceInfo::sys_center)
        .def_readwrite("crystal_dim", &::pygcmc::model::Param::SpaceInfo::crystal_dim)
        .def_readwrite("box_size", &::pygcmc::model::Param::SpaceInfo::box_size)
        .def_readwrite("volume", &::pygcmc::model::Param::SpaceInfo::volume)
        .def_readwrite("target_volume", &::pygcmc::model::Param::SpaceInfo::target_volume)
        .def_readwrite("sys_box_volume", &::pygcmc::model::Param::SpaceInfo::sys_box_volume)
        .def_readwrite("gcmc_volume", &::pygcmc::model::Param::SpaceInfo::gcmc_volume)
        .def_readwrite("protein_volume", &::pygcmc::model::Param::SpaceInfo::protein_volume)
        .def_readwrite("use_vdw_radius_for_grid", &::pygcmc::model::Param::SpaceInfo::use_vdw_radius_for_grid)
        .def_readwrite("exclude_hydrogens_from_grid", &::pygcmc::model::Param::SpaceInfo::exclude_hydrogens_from_grid)
        .def_readwrite("exclude_protein_volume", &::pygcmc::model::Param::SpaceInfo::exclude_protein_volume)
        .def_readwrite("tmp_prob", &::pygcmc::model::Param::SpaceInfo::tmp_prob)
        .def_readwrite("cutoff", &::pygcmc::model::Param::SpaceInfo::cutoff);

    // Bind MCInfo struct
    py::class_<::pygcmc::model::Param::MCInfo>(param, "MCInfo")
        .def(py::init<>())
        .def_readwrite("mc_steps", &::pygcmc::model::Param::MCInfo::mc_steps)
        .def_readwrite("current_step", &::pygcmc::model::Param::MCInfo::current_step)
        .def_readwrite("print_freq", &::pygcmc::model::Param::MCInfo::print_freq)
        .def_readwrite("temperature", &::pygcmc::model::Param::MCInfo::temperature)
        .def_readwrite("beta", &::pygcmc::model::Param::MCInfo::beta)
        .def_readwrite("insertion_deletion_frac", &::pygcmc::model::Param::MCInfo::insertion_deletion_frac)
        .def_readwrite("translation_rotation_frac", &::pygcmc::model::Param::MCInfo::translation_rotation_frac)
        .def_readwrite("max_translation_dist", &::pygcmc::model::Param::MCInfo::max_translation_dist)
        .def_readwrite("max_rotation_angle", &::pygcmc::model::Param::MCInfo::max_rotation_angle)
        .def_readwrite("operation_types", &::pygcmc::model::Param::MCInfo::operation_types)
        .def_readwrite("mc_time_list", &::pygcmc::model::Param::MCInfo::mc_time_list)
        .def_readwrite("mc_time_cumulative", &::pygcmc::model::Param::MCInfo::mc_time_cumulative)
        .def_readwrite("fragment_prob", &::pygcmc::model::Param::MCInfo::fragment_prob)
        .def_readwrite("water_prob", &::pygcmc::model::Param::MCInfo::water_prob)
        .def_readwrite("atom_prob", &::pygcmc::model::Param::MCInfo::atom_prob)
        .def_readwrite("test_prob", &::pygcmc::model::Param::MCInfo::test_prob)
        .def_readwrite("rotate_dih_status", &::pygcmc::model::Param::MCInfo::rotate_dih_status)
        .def_readonly("BOLTZMANN", &::pygcmc::model::Param::MCInfo::BOLTZMANN)
        .def_readonly("KCAL_TO_KJ", &::pygcmc::model::Param::MCInfo::KCAL_TO_KJ);

    // Bind EnergyInfo struct
    py::class_<::pygcmc::model::Param::EnergyInfo>(param, "EnergyInfo")
        .def(py::init<>())
        .def_readwrite("use_group_cutoff", &::pygcmc::model::Param::EnergyInfo::use_group_cutoff)
        .def_readwrite("fragment_cutoff", &::pygcmc::model::Param::EnergyInfo::fragment_cutoff)
        .def_readwrite("protein_cutoff", &::pygcmc::model::Param::EnergyInfo::protein_cutoff)
        .def_readwrite("fragment_cutoff_squared", &::pygcmc::model::Param::EnergyInfo::fragment_cutoff_squared)
        .def_readwrite("protein_cutoff_squared", &::pygcmc::model::Param::EnergyInfo::protein_cutoff_squared)
        .def_readwrite("pairlist_cutoff", &::pygcmc::model::Param::EnergyInfo::pairlist_cutoff)
        .def_readwrite("pairlist_cutoff_squared", &::pygcmc::model::Param::EnergyInfo::pairlist_cutoff_squared)
        .def_readwrite("pairlist_freq", &::pygcmc::model::Param::EnergyInfo::pairlist_freq)
        .def_readwrite("use_switching", &::pygcmc::model::Param::EnergyInfo::use_switching)
        .def_readwrite("switch_dist_fragment", &::pygcmc::model::Param::EnergyInfo::switch_dist_fragment)
        .def_readwrite("switch_dist_protein", &::pygcmc::model::Param::EnergyInfo::switch_dist_protein)
        .def_readwrite("switch_dist_fragment_squared", &::pygcmc::model::Param::EnergyInfo::switch_dist_fragment_squared)
        .def_readwrite("switch_dist_protein_squared", &::pygcmc::model::Param::EnergyInfo::switch_dist_protein_squared)
        .def_readwrite("energy_sw_ref", &::pygcmc::model::Param::EnergyInfo::energy_sw_ref)
        .def_readwrite("energy_sw_scale", &::pygcmc::model::Param::EnergyInfo::energy_sw_scale)
        .def_readwrite("test_sw_filters", &::pygcmc::model::Param::EnergyInfo::test_sw_filters)
        .def_readwrite("apply_sw_filters", &::pygcmc::model::Param::EnergyInfo::apply_sw_filters)
        .def_readwrite("test_energy", &::pygcmc::model::Param::EnergyInfo::test_energy)
        .def_readwrite("pair_list_cutoff_fragment", &::pygcmc::model::Param::EnergyInfo::pair_list_cutoff_fragment)
        .def_readwrite("pair_list_cutoff_protein", &::pygcmc::model::Param::EnergyInfo::pair_list_cutoff_protein)
        .def_readwrite("pair_list_cutoff_fragment_squared", &::pygcmc::model::Param::EnergyInfo::pair_list_cutoff_fragment_squared)
        .def_readwrite("pair_list_cutoff_protein_squared", &::pygcmc::model::Param::EnergyInfo::pair_list_cutoff_protein_squared);

    // Bind FragmentInfo struct
    py::class_<::pygcmc::model::Param::FragmentInfo>(param, "FragmentInfo")
        .def(py::init<>())
        .def_readwrite("water_density", &::pygcmc::model::Param::FragmentInfo::water_density)
        .def_readwrite("epsilon", &::pygcmc::model::Param::FragmentInfo::epsilon)
        .def_readwrite("num_waters", &::pygcmc::model::Param::FragmentInfo::num_waters)
        .def_readwrite("target_num_waters", &::pygcmc::model::Param::FragmentInfo::target_num_waters)
        .def_readwrite("water_index", &::pygcmc::model::Param::FragmentInfo::water_index)
        .def_readwrite("excess_threshold", &::pygcmc::model::Param::FragmentInfo::excess_threshold)
        .def_readwrite("use_number_water_nbar", &::pygcmc::model::Param::FragmentInfo::use_number_water_nbar)
        .def_readwrite("use_const_water_nbar", &::pygcmc::model::Param::FragmentInfo::use_const_water_nbar)
        .def_readwrite("const_water_nbar", &::pygcmc::model::Param::FragmentInfo::const_water_nbar)
        .def_readwrite("init_cutoff", &::pygcmc::model::Param::FragmentInfo::init_cutoff)
        .def_readwrite("init_cutoff_squared", &::pygcmc::model::Param::FragmentInfo::init_cutoff_squared)
        .def_readwrite("use_gcmc_cutoff", &::pygcmc::model::Param::FragmentInfo::use_gcmc_cutoff)
        .def_readwrite("gcmc_cutoff", &::pygcmc::model::Param::FragmentInfo::gcmc_cutoff)
        .def_readwrite("gcmc_cutoff_squared", &::pygcmc::model::Param::FragmentInfo::gcmc_cutoff_squared)
        .def_readwrite("remove_init", &::pygcmc::model::Param::FragmentInfo::remove_init)
        .def_readwrite("remove_excess", &::pygcmc::model::Param::FragmentInfo::remove_excess)
        .def_readwrite("confs_list", &::pygcmc::model::Param::FragmentInfo::confs_list)
        .def_readwrite("cavity_index_list", &::pygcmc::model::Param::FragmentInfo::cavity_index_list)
        .def_readwrite("cavity_list", &::pygcmc::model::Param::FragmentInfo::cavity_list)
        .def_readwrite("cavity_grid_dx_list", &::pygcmc::model::Param::FragmentInfo::cavity_grid_dx_list)
        .def_readwrite("cavity_probe_radius_list", &::pygcmc::model::Param::FragmentInfo::cavity_probe_radius_list)
        .def_readwrite("cavity_mask_list", &::pygcmc::model::Param::FragmentInfo::cavity_mask_list)
        .def_readwrite("conc_list", &::pygcmc::model::Param::FragmentInfo::conc_list)
        .def_readwrite("muex_list", &::pygcmc::model::Param::FragmentInfo::muex_list)
        .def_readwrite("radius_list", &::pygcmc::model::Param::FragmentInfo::radius_list)
        .def_readwrite("conf_list", &::pygcmc::model::Param::FragmentInfo::conf_list)
        .def_readwrite("flag_remove_init", &::pygcmc::model::Param::FragmentInfo::flag_remove_init)
        .def_readwrite("flag_remove_excess", &::pygcmc::model::Param::FragmentInfo::flag_remove_excess)
        .def_readwrite("total_protitp_size", &::pygcmc::model::Param::FragmentInfo::total_protitp_size)
        .def_readwrite("fragconf_list", &::pygcmc::model::Param::FragmentInfo::fragconf_list);

    // Bind BiasInfo struct
    py::class_<::pygcmc::model::Param::BiasInfo>(param, "BiasInfo")
        .def(py::init<>())
        .def_readwrite("use_cavity_bias", &::pygcmc::model::Param::BiasInfo::use_cavity_bias)
        .def_readwrite("sigma", &::pygcmc::model::Param::BiasInfo::sigma)
        .def_readwrite("sigma_squared", &::pygcmc::model::Param::BiasInfo::sigma_squared)
        .def_readwrite("use_conf_bias", &::pygcmc::model::Param::BiasInfo::use_conf_bias)
        .def_readwrite("num_conf_bias_trials", &::pygcmc::model::Param::BiasInfo::num_conf_bias_trials);

    // Bind FileInfo struct
    py::class_<::pygcmc::model::Param::FileInfo>(param, "FileInfo")
        .def(py::init<>())
        .def_readwrite("topology_file", &::pygcmc::model::Param::FileInfo::topology_file)
        .def_readwrite("input_pdb_file", &::pygcmc::model::Param::FileInfo::input_pdb_file)
        .def_readwrite("output_pdb_file", &::pygcmc::model::Param::FileInfo::output_pdb_file)
        .def_readwrite("output_top_file", &::pygcmc::model::Param::FileInfo::output_top_file)
        .def_readwrite("atomtype_file", &::pygcmc::model::Param::FileInfo::atomtype_file)
        .def_readwrite("monomer_dir", &::pygcmc::model::Param::FileInfo::monomer_dir)
        .def_readwrite("conc_norm", &::pygcmc::model::Param::FileInfo::conc_norm)
        .def_readwrite("conc_region", &::pygcmc::model::Param::FileInfo::conc_region)
        .def_readwrite("par_files", &::pygcmc::model::Param::FileInfo::par_files)
        .def_readwrite("protein_top_files", &::pygcmc::model::Param::FileInfo::protein_top_files)
        .def_readwrite("fragment_top_files", &::pygcmc::model::Param::FileInfo::fragment_top_files)
        .def_readwrite("fragment_names", &::pygcmc::model::Param::FileInfo::fragment_names)
        .def_readwrite("fragment_mqtr_files", &::pygcmc::model::Param::FileInfo::fragment_mqtr_files)
        .def_readwrite("tmp_frag_name", &::pygcmc::model::Param::FileInfo::tmp_frag_name)
        .def_readwrite("generate_maps", &::pygcmc::model::Param::FileInfo::generate_maps)
        .def_readwrite("map_prefix", &::pygcmc::model::Param::FileInfo::map_prefix);
}

} // namespace model
} // namespace bindings
} // namespace pygcmc
