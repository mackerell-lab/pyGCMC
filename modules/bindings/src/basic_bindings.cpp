#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/io/ff_parser.hpp"
#include "pygcmc/core/project_atom.hpp"
#include "pygcmc/core/project_residue.hpp"
#include "pygcmc/core/system.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

void init_basic_bindings(py::module& m) {
    // Bind ForceFieldPair struct
    py::class_<io::ForceFieldPair>(m, "ForceFieldPair")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def_readwrite("rmin", &io::ForceFieldPair::rmin)
        .def_readwrite("epsilon", &io::ForceFieldPair::epsilon);

    // Bind PDBAtom struct (for parser usage)
    py::class_<io::PDBAtom>(m, "PDBAtom")
        .def(py::init<>())
        .def_readwrite("serial", &io::PDBAtom::serial)
        .def_readwrite("name", &io::PDBAtom::name)
        .def_readwrite("residue", &io::PDBAtom::residue)
        .def_readwrite("sequence", &io::PDBAtom::sequence)
        .def_readwrite("chain", &io::PDBAtom::chain)
        .def_readwrite("alt_loc", &io::PDBAtom::alt_loc)
        .def_readwrite("insertion_code", &io::PDBAtom::insertion_code)
        .def_readwrite("x", &io::PDBAtom::x)
        .def_readwrite("y", &io::PDBAtom::y)
        .def_readwrite("z", &io::PDBAtom::z)
        .def_readwrite("occupancy", &io::PDBAtom::occupancy)
        .def_readwrite("temp_factor", &io::PDBAtom::temp_factor)
        .def_readwrite("element", &io::PDBAtom::element)
        .def_readwrite("charge", &io::PDBAtom::charge)
        .def_readwrite("type", &io::PDBAtom::type)
        .def_readwrite("topo_type", &io::PDBAtom::topo_type)
        .def_readwrite("topo_charge", &io::PDBAtom::topo_charge)
        .def_readwrite("topo_mass", &io::PDBAtom::topo_mass)
        .def_readwrite("forcefield_epsilon", &io::PDBAtom::forcefield_epsilon)
        .def_readwrite("forcefield_rmin", &io::PDBAtom::forcefield_rmin)
        .def("is_valid", &io::PDBAtom::is_valid)
        .def("has_topology_info", [](const io::PDBAtom& atom) {
            return !atom.topo_type.empty() && !std::isnan(atom.topo_charge) && !std::isnan(atom.topo_mass);
        })
        .def("has_forcefield_info", [](const io::PDBAtom& atom) {
            return !std::isnan(atom.forcefield_epsilon) && !std::isnan(atom.forcefield_rmin);
        });

    // Bind ProjectAtom class
    py::class_<ProjectAtom>(m, "ProjectAtom")
        .def(py::init<>())
        .def(py::init<const io::PDBAtom&>())
        .def_property("serial", &ProjectAtom::get_serial, &ProjectAtom::set_serial)
        .def_property("name", &ProjectAtom::get_name, &ProjectAtom::set_name)
        .def_property("residue", &ProjectAtom::get_residue, &ProjectAtom::set_residue)
        .def_property("sequence", &ProjectAtom::get_sequence, &ProjectAtom::set_sequence)
        .def_property("chain", &ProjectAtom::get_chain, &ProjectAtom::set_chain)
        .def_property("alt_loc", &ProjectAtom::get_alt_loc, &ProjectAtom::set_alt_loc)
        .def_property("insertion_code", &ProjectAtom::get_insertion_code, &ProjectAtom::set_insertion_code)
        .def_property("x", &ProjectAtom::get_x, &ProjectAtom::set_x)
        .def_property("y", &ProjectAtom::get_y, &ProjectAtom::set_y)
        .def_property("z", &ProjectAtom::get_z, &ProjectAtom::set_z)
        .def_property("occupancy", &ProjectAtom::get_occupancy, &ProjectAtom::set_occupancy)
        .def_property("temp_factor", &ProjectAtom::get_temp_factor, &ProjectAtom::set_temp_factor)
        .def_property("element", &ProjectAtom::get_element, &ProjectAtom::set_element)
        .def_property("charge", &ProjectAtom::get_charge, &ProjectAtom::set_charge)
        .def_property("type", &ProjectAtom::get_type, &ProjectAtom::set_type)
        .def_property("topo_type", &ProjectAtom::get_topo_type, &ProjectAtom::set_topo_type)
        .def_property("topo_charge", &ProjectAtom::get_topo_charge, &ProjectAtom::set_topo_charge)
        .def_property("topo_mass", &ProjectAtom::get_topo_mass, &ProjectAtom::set_topo_mass)
        .def_property("forcefield_epsilon", &ProjectAtom::get_forcefield_epsilon, &ProjectAtom::set_forcefield_epsilon)
        .def_property("forcefield_rmin", &ProjectAtom::get_forcefield_rmin, &ProjectAtom::set_forcefield_rmin)
        .def("is_valid", &ProjectAtom::is_valid)
        .def("has_topology_info", &ProjectAtom::has_topology_info)
        .def("has_forcefield_info", &ProjectAtom::has_forcefield_info);

    // Bind ProjectResidue class
    py::class_<ProjectResidue>(m, "ProjectResidue")
        .def(py::init<>())
        .def(py::init<const std::string&, int, char>())
        .def_property("name", &ProjectResidue::get_name, &ProjectResidue::set_name)
        .def_property("sequence_number", &ProjectResidue::get_sequence_number, &ProjectResidue::set_sequence_number)
        .def_property("chain_id", &ProjectResidue::get_chain_id, &ProjectResidue::set_chain_id)
        .def_property("atoms", &ProjectResidue::get_atoms, &ProjectResidue::set_atoms)
        .def("center_of_mass", &ProjectResidue::center_of_mass)
        .def("atom_count", &ProjectResidue::atom_count);

    // Bind IOResidue struct with shared_ptr (for parser usage)
    py::class_<io::IOResidue, std::shared_ptr<io::IOResidue>>(m, "IOResidue")
        .def(py::init<>())
        .def_readwrite("name", &io::IOResidue::name)
        .def_readwrite("sequence_number", &io::IOResidue::sequence_number)
        .def_readwrite("chain_id", &io::IOResidue::chain_id)
        .def_readwrite("atoms", &io::IOResidue::atoms)
        .def_readwrite("atom_ptrs", &io::IOResidue::atom_ptrs)
        .def("center_of_mass", &io::IOResidue::center_of_mass)
        .def("atom_count", &io::IOResidue::atom_count);

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
        .def("get_particle_count", &Residue::get_particle_count)
        .def("center_of_mass", &Residue::center_of_mass);
} 