#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pygcmc/core/system.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

void init_system_bindings(py::module& m) {
    // Bind System class
    py::class_<System>(m, "System")
        .def(py::init<const std::string&, const std::string&>(),
             py::kw_only(),
             py::arg("pdb") = "", py::arg("top") = "",
             "Create a system, optionally loading structure from PDB and topology files")
        .def("load_structure", &System::load_structure,
             py::arg("pdb_file"), py::arg("top_file"),
             "Load structure from PDB and topology files")
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