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
#include "pygcmc/core/project_atom.hpp"
#include "pygcmc/core/project_residue.hpp"

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

    // Bind PDBParser class
    py::class_<io::PDBParser>(m, "PDBParser")
        .def(py::init<>())
        .def_static("parse", &io::PDBParser::parse);

    // Bind Residue struct for PDB with proper atom handling
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
        .def("update_pdb_atoms", [](io::TopParser& self, py::list atoms) -> int {
            std::vector<io::PDBAtom*> c_atoms;
            c_atoms.reserve(py::len(atoms));
            
            for(auto item : atoms) {
                try {
                    // Try as raw PDBAtom first
                    auto& atom = item.cast<io::PDBAtom&>();
                    c_atoms.push_back(&atom);
                } catch (const py::cast_error&) {
                    try {
                        // Try as ProjectAtom
                        auto& wrapper = item.cast<ProjectAtom&>();
                        c_atoms.push_back(wrapper.get_ptr().get());
                    } catch (const py::cast_error&) {
                        // Try as shared_ptr
                        auto shared_atom = item.cast<std::shared_ptr<io::PDBAtom>>();
                        c_atoms.push_back(shared_atom.get());
                    }
                }
            }
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
            atom_ptrs.reserve(py::len(atoms));
            
            for (auto item : atoms) {
                try {
                    // Try as raw PDBAtom first
                    auto& atom = item.cast<io::PDBAtom&>();
                    atom_ptrs.push_back(&atom);
                } catch (const py::cast_error&) {
                    try {
                        // Try as ProjectAtom
                        auto& wrapper = item.cast<ProjectAtom&>();
                        atom_ptrs.push_back(wrapper.get_ptr().get());
                    } catch (const py::cast_error&) {
                        // Try as shared_ptr
                        auto shared_atom = item.cast<std::shared_ptr<io::PDBAtom>>();
                        atom_ptrs.push_back(shared_atom.get());
                    }
                }
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
        .def("load_structure", &Project::load_structure, "Load structure from PDB and topology files")
        .def("load_forcefield", &Project::load_forcefield, "Load force field from parameter files")
        .def("get_name", &Project::get_name, "Get project name")
        .def("print_atom_info", &Project::print_atom_info, "Print detailed information for a single atom")
        .def("print_detailed_atom_info", &Project::print_detailed_atom_info, "Print detailed information for a single atom in table format")
        .def("print_atom_table_header", &Project::print_atom_table_header, "Print header for atom table")
        .def("print_all_atoms", &Project::print_all_atoms, "Print information for all atoms")
        .def("print_forcefield_info", &Project::print_forcefield_info, "Print force field information")
        .def("print_nbfix_info", &Project::print_nbfix_info, "Print NBFIX parameters")
        .def("print_global_parameters", &Project::print_global_parameters, "Print global force field parameters");

    // Bind Structure class with proper wrapper support
    py::class_<Structure, std::shared_ptr<Structure>>(m, "Structure")
        .def(py::init<>())
        .def("apply_forcefield", &Structure::apply_forcefield)
        .def_property_readonly("residues", [](Structure& self) {
            std::vector<ProjectResidue> wrapped_residues;
            for (const auto& residue : self.residues()) {
                ProjectResidue wrapper(residue->name, residue->sequence_number, residue->chain_id);
                std::vector<ProjectAtom> atoms;
                for (const auto& atom_ptr : residue->atom_ptrs) {
                    if (atom_ptr) {
                        atoms.emplace_back(*atom_ptr);
                    }
                }
                wrapper.set_atoms(atoms);
                wrapped_residues.push_back(wrapper);
            }
            return wrapped_residues;
        })
        .def_property_readonly("atoms", [](Structure& self) {
            std::vector<ProjectAtom> wrapped_atoms;
            for (const auto& atom_ptr : self.atoms()) {
                if (atom_ptr != nullptr) {
                    wrapped_atoms.emplace_back(*atom_ptr);
                }
            }
            return wrapped_atoms;
        })
        .def("__len__", &Structure::get_num_atoms);

    // Bind ForceField class
    py::class_<ForceField, std::shared_ptr<ForceField>>(m, "ForceField")
        .def(py::init<>())
        .def_property("cutoff", &ForceField::get_cutoff, &ForceField::set_cutoff)
        .def_property("switching", &ForceField::get_switching, &ForceField::set_switching)
        .def_property("pairlist_distance", &ForceField::get_pairlist_distance, &ForceField::set_pairlist_distance)
        .def_property_readonly("nonbonded_params", 
            static_cast<const std::map<std::string, io::ForceFieldPair>& (ForceField::*)() const>(&ForceField::nonbonded_params))
        .def_property_readonly("nbfix_params", 
            static_cast<const std::map<std::pair<std::string, std::string>, io::ForceFieldPair>& (ForceField::*)() const>(&ForceField::nbfix_params))
        .def("print_nonbonded_params", &ForceField::print_nonbonded_params, "Print all nonbonded parameters");
}
