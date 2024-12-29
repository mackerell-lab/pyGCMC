// modules/bindings/src/bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/system.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

PYBIND11_MODULE(pyGCMC_bindings, m) {
    m.doc() = "Python bindings for pyGCMC simulation library";

    // Bind Particle struct
    py::class_<Particle>(m, "Particle")
        .def(py::init<const std::array<double, 3>&, const std::array<double, 3>&, double, double>(),
             py::arg("position") = std::array<double, 3>{0, 0, 0},
             py::arg("velocity") = std::array<double, 3>{0, 0, 0},
             py::arg("charge") = 0.0,
             py::arg("mass") = 1.0)
        .def_readwrite("position", &Particle::position)
        .def_readwrite("velocity", &Particle::velocity)
        .def_readwrite("charge", &Particle::charge)
        .def_readwrite("mass", &Particle::mass)
        .def_readwrite("is_virtual", &Particle::is_virtual)
        .def("is_valid", &Particle::is_valid);

    // Bind Residue struct
    py::class_<Residue>(m, "Residue")
        .def(py::init<const std::string&>(),
             py::arg("name") = "")
        .def_readwrite("name", &Residue::name)
        .def_readwrite("particles", &Residue::particles)
        .def("center_of_mass", &Residue::center_of_mass)
        .def("atom_count", &Residue::atom_count);

    // Bind System class
    py::class_<System>(m, "System")
        .def(py::init<>())
        .def("add_residue", &System::add_residue, py::arg("name"), "Add a residue to the system")
        .def("remove_residue", &System::remove_residue, py::arg("index"), "Remove a residue by index")
        .def("get_residue_count", &System::get_residue_count, "Get number of residues")
        .def("get_residue", 
             py::overload_cast<size_t>(&System::get_residue, py::const_),
             py::arg("index"),
             "Get residue by index (const)")
        .def("get_residue",
             py::overload_cast<size_t>(&System::get_residue),
             py::arg("index"),
             "Get residue by index (mutable)")
        .def("add_particle", &System::add_particle, py::arg("residue_index"), py::arg("particle"),
             "Add a particle to a residue")
        .def("remove_particle", &System::remove_particle, py::arg("residue_index"), py::arg("particle_index"),
             "Remove a particle from a residue")
        .def("get_particle_count", &System::get_particle_count, py::arg("residue_index"),
             "Get number of particles in a residue")
        .def("get_particle",
             py::overload_cast<size_t, size_t>(&System::get_particle, py::const_),
             py::arg("residue_index"), py::arg("particle_index"),
             "Get particle by indices (const)")
        .def("get_particle",
             py::overload_cast<size_t, size_t>(&System::get_particle),
             py::arg("residue_index"), py::arg("particle_index"),
             "Get particle by indices (mutable)")
        .def("get_particle_mass", &System::get_particle_mass,
             py::arg("residue_index"), py::arg("particle_index"),
             "Get mass of a particle")
        .def("set_particle_mass", &System::set_particle_mass,
             py::arg("residue_index"), py::arg("particle_index"), py::arg("mass"),
             "Set mass of a particle")
        .def("set_virtual_site", &System::set_virtual_site,
             py::arg("residue_index"), py::arg("particle_index"), py::arg("is_virtual"),
             "Set whether a particle is a virtual site")
        .def("is_virtual_site", &System::is_virtual_site,
             py::arg("residue_index"), py::arg("particle_index"),
             "Check if a particle is a virtual site")
        .def("add_constraint", &System::add_constraint,
             py::arg("residue1"), py::arg("particle1"), py::arg("residue2"), py::arg("particle2"), py::arg("distance"),
             "Add a constraint between two particles")
        .def("remove_constraint", &System::remove_constraint, py::arg("index"),
             "Remove a constraint by index")
        .def("get_constraint_count", &System::get_constraint_count,
             "Get number of constraints")
        .def("get_constraint", &System::get_constraint, py::arg("index"),
             "Get constraint by index")
        .def("compute_energy", &System::compute_energy,
             "Compute total system energy")
        .def("update_positions", &System::update_positions, py::arg("dt"),
             "Update particle positions")
        .def("update_velocities", &System::update_velocities, py::arg("dt"),
             "Update particle velocities")
        .def("set_periodic_box_vectors", &System::set_periodic_box_vectors,
             py::arg("a"), py::arg("b"), py::arg("c"),
             "Set periodic boundary conditions box vectors")
        .def("get_periodic_box_vectors", &System::get_periodic_box_vectors,
             py::arg("a"), py::arg("b"), py::arg("c"),
             "Get periodic boundary conditions box vectors")
        .def("uses_periodic_boundary_conditions", &System::uses_periodic_boundary_conditions,
             "Check if system uses periodic boundary conditions")
        .def("compute_distance", &System::compute_distance,
             py::arg("p1"), py::arg("p2"),
             "Compute distance between two particles");
}
