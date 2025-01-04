#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pygcmc/core/system.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

void init_system_bindings(py::module& m) {
    // Bind System class
    py::class_<System>(m, "System")
        .def(py::init<>())
        .def(py::init<const std::string&, const std::string&>(),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Create a system and load structure from PDB and PSF files")
        .def(py::init([](const std::string& pdb, const std::vector<std::string>& psf) {
                System system;
                for (const auto& psf_file : psf) {
                    system.load_structure_psf(pdb, psf_file);
                }
                return system;
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Create a system and load structure from PDB and multiple PSF files")
        .def(py::init([](const std::string& pdb, const std::string& top) {
                return System::from_top(pdb, top);
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("top"),
             "Create a system and load structure from PDB and TOP files")
        .def_static("from_top", &System::from_top,
             py::kw_only(),
             py::arg("pdb"), py::arg("top"),
             "Create a system and load structure from PDB and TOP files")
        .def("load_structure",
             [](System& self, const std::string& pdb, const std::string& psf) {
                 self.load_structure_psf(pdb, psf);
             },
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Load structure from PDB and PSF files")
        .def("load_structure",
             [](System& self, const std::string& pdb, const std::vector<std::string>& psf) {
                 for (const auto& psf_file : psf) {
                     self.load_structure_psf(pdb, psf_file);
                 }
             },
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Load structure from PDB and multiple PSF files")
        .def("load_structure",
             [](System& self, const std::string& pdb, const std::string& top) {
                 self.load_structure_top(pdb, top);
             },
             py::kw_only(),
             py::arg("pdb"), py::arg("top"),
             "Load structure from PDB and TOP files")
        .def("load_structure",
             static_cast<void (System::*)(const Structure&)>(&System::load_structure),
             py::arg("structure"),
             "Load structure from a Structure object")
        .def("load_structure_psf_auto",
             &System::load_structure_psf_auto,
             py::arg("pdb_file"), py::arg("psf_file"),
             "Load structure from PDB and PSF files with automatic detection of PSF type")
        .def("load_structure_psf_multi",
             &System::load_structure_psf_multi,
             py::arg("pdb_file"), py::arg("psf_file"),
             "Load structure from PDB and multi-residue PSF file")
        .def("load_structure_psf_single",
             &System::load_structure_psf_single,
             py::arg("pdb_file"), py::arg("psf_file"), py::arg("target_residue"),
             "Load structure from PDB and single-residue PSF file, applying to specified residue type")
        // Residue management
        .def("get_residue_count", &System::get_residue_count,
             "Get number of residues")
        .def("add_residue", &System::add_residue,
             py::arg("name"),
             "Add a residue to the system")
        .def("remove_residue", &System::remove_residue,
             py::arg("index"),
             "Remove a residue by index")
        .def("get_residue",
             py::overload_cast<size_t>(&System::get_residue, py::const_),
             py::arg("index"),
             "Get residue by index (const)")
        .def("get_residue",
             py::overload_cast<size_t>(&System::get_residue),
             py::arg("index"),
             "Get residue by index (mutable)")
        // Particle management
        .def("get_particle_count", &System::get_particle_count,
             py::arg("residue_index"),
             "Get number of particles in a residue")
        .def("add_particle", &System::add_particle,
             py::arg("residue_index"), py::arg("particle"),
             "Add a particle to a residue")
        .def("remove_particle", &System::remove_particle,
             py::arg("residue_index"), py::arg("particle_index"),
             "Remove a particle from a residue")
        .def("get_particle",
             py::overload_cast<size_t, size_t>(&System::get_particle, py::const_),
             py::arg("residue_index"), py::arg("particle_index"),
             "Get particle by indices (const)")
        .def("get_particle",
             py::overload_cast<size_t, size_t>(&System::get_particle),
             py::arg("residue_index"), py::arg("particle_index"),
             "Get particle by indices (mutable)")
        // PDBAtom support
        .def("get_pdb_atom_count", &System::get_pdb_atom_count,
             "Get number of PDB atoms")
        .def("get_pdb_atom",
             py::overload_cast<size_t>(&System::get_pdb_atom, py::const_),
             py::arg("index"),
             "Get PDB atom by index (const)")
        .def("get_pdb_atom",
             py::overload_cast<size_t>(&System::get_pdb_atom),
             py::arg("index"),
             "Get PDB atom by index (mutable)")
        .def("add_pdb_atom", &System::add_pdb_atom,
             py::arg("atom"),
             "Add a PDB atom to the system")
        .def("remove_pdb_atom", &System::remove_pdb_atom,
             py::arg("index"),
             "Remove a PDB atom by index")
        .def("get_pdb_atoms_by_residue", &System::get_pdb_atoms_by_residue,
             py::arg("residue_name"),
             "Get PDB atoms by residue name")
        .def("get_pdb_atoms_by_residue_sequence", &System::get_pdb_atoms_by_residue_sequence,
             py::arg("residue_name"), py::arg("sequence"),
             "Get PDB atoms by residue name and sequence number")
        .def("get_pdb_atoms_by_chain", &System::get_pdb_atoms_by_chain,
             py::arg("chain"),
             "Get PDB atoms by chain identifier")
        .def("clear_pdb_atoms", &System::clear_pdb_atoms,
             "Clear all PDB atoms")
        .def("has_pdb_atoms", &System::has_pdb_atoms,
             "Check if system has any PDB atoms");
} 