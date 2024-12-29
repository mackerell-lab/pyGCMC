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
        .def(py::init<int, const std::string&, const std::string&, 
                 int, double, double, double, 
                 double, const std::string&, const std::string&>(),
             py::arg("serial") = 0,
             py::arg("name") = "",
             py::arg("residue") = "",
             py::arg("sequence") = 0,
             py::arg("x") = 0.0,
             py::arg("y") = 0.0,
             py::arg("z") = 0.0,
             py::arg("charge") = 0.0,
             py::arg("type") = "",
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
        .def("position", &Particle::position, "Get particle position as [x, y, z] array")
        .def("velocity", &Particle::velocity, "Get particle velocity as [vx, vy, vz] array")
        .def("set_position", &Particle::set_position, py::arg("pos"), 
             "Set particle position from [x, y, z] array")
        .def("set_velocity", &Particle::set_velocity, py::arg("vel"), 
             "Set particle velocity from [vx, vy, vz] array");

    // Bind Residue struct
    py::class_<Residue>(m, "Residue")
        .def(py::init<const std::string&, int, char>(),
             py::arg("name") = "",
             py::arg("sequence_number") = 0,
             py::arg("chain_id") = ' ')
        .def_readwrite("name", &Residue::name)
        .def_readwrite("sequence_number", &Residue::sequence_number)
        .def_readwrite("chain_id", &Residue::chain_id)
        .def_readwrite("atoms", &Residue::atoms)
        .def("center_of_mass", &Residue::center_of_mass, 
             "Compute geometric center of the residue (assuming equal masses)")
        .def("atom_count", &Residue::atom_count, 
             "Get the number of atoms in this residue");

    // Bind System class
    py::class_<System>(m, "System")
        .def(py::init<double, double>(),
             py::arg("epsilon") = 1.0,
             py::arg("sigma") = 1.0)
        .def("add_residue", &System::add_residue, py::arg("residue"), "Add a residue to the system")
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
        .def("compute_total_energy", &System::compute_total_energy, "Compute total system energy")
        .def("get_system_state", &System::get_system_state, "Get system kinetic and potential energy")
        .def("update_positions", &System::update_positions, py::arg("dt"), "Update particle positions")
        .def("update_velocities", &System::update_velocities, py::arg("dt"), "Update particle velocities")
        .def("set_periodic_boundary", &System::set_periodic_boundary, py::arg("box_size"), "Set periodic boundary conditions")
        .def("load_pdb", &System::load_pdb, py::arg("filename"), "Load particles from PDB file")
        .def("load_psf", &System::load_psf, py::arg("filename"), "Load system topology from PSF file")
        .def("load_top", &System::load_top, py::arg("filename"), "Load force field parameters from TOP file")
        .def("load_itp", &System::load_itp, py::arg("filename"), "Load molecular information from ITP file")
        .def("load_forcefield", &System::load_forcefield, py::arg("filename"), "Load non-bonded parameters from force field file");
}
