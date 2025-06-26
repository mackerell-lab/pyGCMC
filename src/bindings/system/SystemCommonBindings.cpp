// src/bindings/system/SystemCommonBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "system/SystemModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace system {

void init_common_bindings(py::module& m, py::module&) {
    // Add LogLevel enum
    py::enum_<::pygcmc::system::common::LogLevel>(m, "LogLevel")
        .value("DEBUG", ::pygcmc::system::common::LogLevel::DEBUG)
        .value("INFO", ::pygcmc::system::common::LogLevel::INFO)
        .value("WARNING", ::pygcmc::system::common::LogLevel::WARNING)
        .value("ERROR", ::pygcmc::system::common::LogLevel::ERROR);

    py::class_<::pygcmc::system::System>(m, "System")
        .def(py::init<>())
        .def_static("set_verbose", &::pygcmc::system::System::set_verbose, 
                   "Set verbose mode for logging")
        .def_static("set_log_level", &::pygcmc::system::System::set_log_level, 
                   "Set the minimum log level")
        .def("initialize_parameters", &::pygcmc::system::System::initialize_parameters, 
             "Initialize all system parameters")
        .def("process_cavity_list", &::pygcmc::system::System::process_cavity_list,
             "Process and sort the cavity list")
        .def("initialize_mc_time_list", &::pygcmc::system::System::initialize_mc_time_list,
             "Initialize Monte Carlo time list");

    // Bind TypeMaps
    py::class_<pygcmc::model::TypeMaps>(m, "TypeMaps")
        .def(py::init<>())
        .def("get_or_add_type", &pygcmc::model::TypeMaps::getOrAddType)
        .def("get_type_name", &pygcmc::model::TypeMaps::getTypeName)
        .def_readonly("atomTypes", &pygcmc::model::TypeMaps::atomTypes);
}

} // namespace system
} // namespace bindings
} // namespace pygcmc