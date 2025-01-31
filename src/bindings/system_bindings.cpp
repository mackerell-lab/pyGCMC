#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "system/system.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void bind_system(py::module& m) {
    py::class_<system::System>(m, "System")
        .def(py::init<>())
        .def("initialize_parameters", &system::System::initialize_parameters, 
             "Initialize all system parameters")
        .def("process_cavity_list", &system::System::process_cavity_list,
             "Process and sort the cavity list")
        .def("initialize_mc_time_list", &system::System::initialize_mc_time_list,
             "Initialize Monte Carlo time list");
}

} // namespace bindings
} // namespace pygcmc
