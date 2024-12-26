// modules/core/src/bindings.cpp
#include <pybind11/pybind11.h>
#include "pygcmc/core/system.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyGCMC, m) {
    py::class_<Particle>(m, "Particle")
        .def(py::init<double, double, double, int, double>(),
             py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0,
             py::arg("type") = 0, py::arg("charge") = 0.0)
        .def_readwrite("x", &Particle::x)
        .def_readwrite("y", &Particle::y)
        .def_readwrite("z", &Particle::z)
        .def_readwrite("type", &Particle::type)
        .def_readwrite("charge", &Particle::charge);

    py::class_<System>(m, "System")
        .def(py::init<double, double>(), py::arg("epsilon") = 1.0, py::arg("sigma") = 1.0)
        .def("add_particle", &System::add_particle)
        .def("remove_particle", &System::remove_particle)
        .def("compute_total_energy", &System::compute_total_energy)
        .def("get_particle_count", &System::get_particle_count);
}
