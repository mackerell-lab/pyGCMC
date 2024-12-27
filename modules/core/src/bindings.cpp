// modules/core/src/io/bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/system.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

PYBIND11_MODULE(pyGCMC, m) {
    m.doc() = "Python bindings for GCMC simulation library";

    py::class_<Particle>(m, "Particle")
        .def(py::init<int, const std::string&, const std::string&, 
                     int, double, double, double, 
                     double, int, const std::string&>(),
             py::arg("serial") = 0,
             py::arg("name") = "",
             py::arg("residue") = "",
             py::arg("sequence") = 0,
             py::arg("x") = 0.0,
             py::arg("y") = 0.0,
             py::arg("z") = 0.0,
             py::arg("charge") = 0.0,
             py::arg("type") = 0,
             py::arg("nameTop") = "")
        .def_readwrite("serial", &Particle::serial)
        .def_readwrite("name", &Particle::name)
        .def_readwrite("residue", &Particle::residue)
        .def_readwrite("sequence", &Particle::sequence)
        .def_readwrite("x", &Particle::x)
        .def_readwrite("y", &Particle::y)
        .def_readwrite("z", &Particle::z)
        .def_readwrite("charge", &Particle::charge)
        .def_readwrite("type", &Particle::type)
        .def_readwrite("nameTop", &Particle::nameTop)
        .def_readwrite("typeNum", &Particle::typeNum)
        .def_readwrite("vx", &Particle::vx)
        .def_readwrite("vy", &Particle::vy)
        .def_readwrite("vz", &Particle::vz)
        .def("is_valid", &Particle::is_valid)
        .def("position", &Particle::position)
        .def("velocity", &Particle::velocity)
        .def("set_position", &Particle::set_position)
        .def("set_velocity", &Particle::set_velocity);

    py::class_<System>(m, "System")
        .def(py::init<double, double>(),
             py::arg("epsilon") = 1.0,
             py::arg("sigma") = 1.0)
        .def("add_particle", &System::add_particle,
             py::arg("particle"),
             "Add a particle to the system")
        .def("remove_particle", &System::remove_particle,
             py::arg("index"),
             "Remove a particle by index")
        .def("compute_total_energy", &System::compute_total_energy,
             "Compute total system energy")
        .def("get_particle_count", &System::get_particle_count,
             "Get number of particles")
        .def("get_particle", 
             py::overload_cast<size_t>(&System::get_particle, py::const_),
             py::arg("index"),
             "Get particle by index (const)")
        .def("get_particle",
             py::overload_cast<size_t>(&System::get_particle),
             py::arg("index"),
             "Get particle by index (mutable)")
        .def("update_positions", &System::update_positions,
             py::arg("dt"),
             "Update particle positions")
        .def("update_velocities", &System::update_velocities,
             py::arg("dt"),
             "Update particle velocities")
        .def("set_periodic_boundary", &System::set_periodic_boundary,
             py::arg("box_size"),
             "Set periodic boundary conditions")
        .def("get_system_state", &System::get_system_state,
             "Get system kinetic and potential energy")
        .def("load_pdb", &System::load_pdb,
             py::arg("filename"),
             "Load particles from PDB file")
        .def("load_psf", &System::load_psf,
             py::arg("filename"),
             "Load system topology from PSF file")
        .def("load_top", &System::load_top,
             py::arg("filename"),
             "Load force field parameters from TOP file")
        .def("load_itp", &System::load_itp,
             py::arg("filename"),
             "Load molecular information from ITP file")
        .def("load_forcefield", &System::load_forcefield,
             py::arg("filename"),
             "Load non-bonded parameters from force field file");
}
