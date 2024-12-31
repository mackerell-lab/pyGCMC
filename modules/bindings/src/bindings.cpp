// modules/bindings/src/bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pygcmc/core/system.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/ff_parser.hpp"
#include "pygcmc/core/project.hpp"
#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for pygcmc simulation library";

    // Bind ForceFieldPair struct
    py::class_<io::ForceFieldPair>(m, "ForceFieldPair")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def_readwrite("rmin", &io::ForceFieldPair::rmin)
        .def_readwrite("epsilon", &io::ForceFieldPair::epsilon);

    // Bind PDBAtom struct
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

    // Bind IOResidue struct
    py::class_<io::IOResidue>(m, "IOResidue")
        .def(py::init<>())
        .def_readwrite("name", &io::IOResidue::name)
        .def_readwrite("sequence_number", &io::IOResidue::sequence_number)
        .def_readwrite("chain_id", &io::IOResidue::chain_id)
        .def_readwrite("atoms", &io::IOResidue::atoms)
        .def_readwrite("atom_ptrs", &io::IOResidue::atom_ptrs)
        .def("center_of_mass", &io::IOResidue::center_of_mass)
        .def("atom_count", &io::IOResidue::atom_count);

    // Bind PDBParser class
    py::class_<io::PDBParser>(m, "PDBParser")
        .def(py::init<>())
        .def_static("parse", &io::PDBParser::parse);

    // Bind Residue struct for PDB
    py::class_<io::Residue>(m, "PDBResidue")
        .def(py::init<>())
        .def_readwrite("name", &io::Residue::name)
        .def_readwrite("sequence_number", &io::Residue::sequence_number)
        .def_readwrite("chain_id", &io::Residue::chain_id)
        .def_readwrite("atoms", &io::Residue::atoms);

    // Bind TopParser class
    py::class_<io::TopParser>(m, "TopParser")
        .def(py::init<>())
        .def("parse", &io::TopParser::parse)
        .def("parse_with_includes", &io::TopParser::parse_with_includes)
        // 绑定新的 update_pdb_atoms 方法，仅接受指针版本
        .def("update_pdb_atoms", [](io::TopParser& self, py::list atoms) -> int {
            std::vector<io::PDBAtom*> c_atoms;
            for(auto item : atoms){
                // 确保 item 是 PDBAtom 的实例
                io::PDBAtom* atom = item.cast<io::PDBAtom*>();
                c_atoms.push_back(atom);
            }
            // 调用 C++ 的 update_pdb_atoms 方法
            return self.update_pdb_atoms(c_atoms);
        }, py::arg("atoms"));

    // Bind FFParser class
    py::class_<io::FFParser>(m, "FFParser")
        .def(py::init<>())
        .def("parse", &io::FFParser::parse)
        .def("get_nonbonded_params", &io::FFParser::get_nonbonded_params)
        .def("get_nbfix_params", &io::FFParser::get_nbfix_params)
        .def("update_pdb_atoms", [](io::FFParser& self, py::list atoms) -> int {
            std::vector<io::PDBAtom*> atom_ptrs;
            for (auto item : atoms) {
                atom_ptrs.push_back(item.cast<io::PDBAtom*>());
            }
            return self.update_pdb_atoms(atom_ptrs);
        })
        // Add global nonbonded parameter getters
        .def("get_cutnb", &io::FFParser::get_cutnb)
        .def("get_ctofnb", &io::FFParser::get_ctofnb)
        .def("get_ctonnb", &io::FFParser::get_ctonnb)
        .def("get_eps", &io::FFParser::get_eps)
        .def("get_e14fac", &io::FFParser::get_e14fac)
        .def("get_wmin", &io::FFParser::get_wmin)
        // Add static method for merging NBFIX parameters
        .def_static("merge_nbfix_params", [](py::list parsers) {
            std::vector<const io::FFParser*> parser_ptrs;
            for (auto item : parsers) {
                parser_ptrs.push_back(item.cast<io::FFParser*>());
            }
            return io::FFParser::merge_nbfix_params(parser_ptrs);
        });

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

    // Bind Project class
    py::class_<Project>(m, "Project")
        .def(py::init<const std::string&>(), py::arg("name") = "")
        .def("load_structure", &Project::load_structure)
        .def("load_forcefield", &Project::load_forcefield)
        .def("get_name", &Project::get_name);

    // Bind Structure class
    py::class_<Structure, std::shared_ptr<Structure>>(m, "Structure")
        .def(py::init<>())
        .def("apply_forcefield", &Structure::apply_forcefield)
        .def_property_readonly("residues", &Structure::residues)
        .def_property_readonly("atoms", &Structure::atoms)
        .def("__len__", &Structure::get_num_atoms);

    // Bind ForceField class
    py::class_<ForceField, std::shared_ptr<ForceField>>(m, "ForceField")
        .def(py::init<>())
        .def_property("cutoff", &ForceField::get_cutoff, &ForceField::set_cutoff)
        .def_property("switching", &ForceField::get_switching, &ForceField::set_switching)
        .def_property("pairlist_distance", &ForceField::get_pairlist_distance, &ForceField::set_pairlist_distance)
        .def_property_readonly("nonbonded_params", &ForceField::nonbonded_params)
        .def_property_readonly("nbfix_params", &ForceField::nbfix_params);
}
