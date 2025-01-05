#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/system.hpp"
#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/io/psf_parser.hpp"
#include <variant>
#include <unordered_map>

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
                return System::from_pdb_psf(pdb, psf);
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Create a system and load structure from PDB and multiple PSF files")
        .def(py::init([](const std::string& pdb, const std::vector<std::string>& psf, const std::string& itp) {
                return System::from_pdb_psf_itp(pdb, psf, itp);
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"), py::arg("itp"),
             "Create a system and load structure from PDB, multiple PSF files, and an ITP file")
        .def(py::init([](const std::string& pdb, const std::vector<std::string>& psf, const std::vector<std::string>& itp) {
                return System::from_pdb_psf_itps(pdb, psf, itp);
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"), py::arg("itp"),
             "Create a system and load structure from PDB, multiple PSF files, and multiple ITP files")
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
             [](System& self, const std::string& pdb, const std::string& top) {
                 self.load_structure_top(pdb, top);
             },
             py::kw_only(),
             py::arg("pdb"), py::arg("top"),
             "Load structure from PDB and TOP files")
        .def("load_structure",
             [](System& self, py::kwargs kwargs) {
                 std::unordered_map<std::string, std::variant<std::string, std::vector<std::string>>> cpp_kwargs;
                 
                 // Convert Python kwargs to C++ map
                 for (const auto& item : kwargs) {
                     std::string key = py::cast<std::string>(item.first);
                     py::handle value = item.second;
                     
                     if (py::isinstance<py::str>(value)) {
                         cpp_kwargs[key] = py::cast<std::string>(value);
                     } else if (py::isinstance<py::list>(value)) {
                         cpp_kwargs[key] = py::cast<std::vector<std::string>>(value);
                     }
                 }
                 
                 // Instead of creating a new system and copying, directly load into this system
                 if (cpp_kwargs.count("pdb") && cpp_kwargs.count("psf")) {
                     const std::string& pdb = std::get<std::string>(cpp_kwargs["pdb"]);
                     if (std::holds_alternative<std::string>(cpp_kwargs["psf"])) {
                         self.load_structure_psf_auto(pdb, std::get<std::string>(cpp_kwargs["psf"]));
                     } else {
                         const auto& psf_files = std::get<std::vector<std::string>>(cpp_kwargs["psf"]);
                         // Use the same method as the constructor
                         System temp = System::from_pdb_psf(pdb, psf_files);
                         self = std::move(temp);
                     }
                 } else if (cpp_kwargs.count("pdb") && cpp_kwargs.count("top")) {
                     const std::string& pdb = std::get<std::string>(cpp_kwargs["pdb"]);
                     const std::string& top = std::get<std::string>(cpp_kwargs["top"]);
                     self.load_structure_top(pdb, top);
                 }
             },
             "Load structure with flexible file loading options.\n\n"
             "Args:\n"
             "    pdb (str): Path to PDB file\n"
             "    psf (str or List[str], optional): Path(s) to PSF file(s)\n"
             "    top (str, optional): Path to TOP file\n"
             "    itp (str or List[str], optional): Path(s) to ITP file(s)\n\n"
             "Note: You can combine different file types as needed.")
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
        .def(py::init([](py::kwargs kwargs) {
                std::unordered_map<std::string, std::variant<std::string, std::vector<std::string>>> cpp_kwargs;
                
                // Convert Python kwargs to C++ map
                for (const auto& item : kwargs) {
                    std::string key = py::cast<std::string>(item.first);
                    py::handle value = item.second;
                    
                    if (py::isinstance<py::str>(value)) {
                        cpp_kwargs[key] = py::cast<std::string>(value);
                    } else if (py::isinstance<py::list>(value)) {
                        cpp_kwargs[key] = py::cast<std::vector<std::string>>(value);
                    }
                }
                
                return System::from_kwargs(cpp_kwargs);
             }),
             "Create a system with flexible file loading options.\n\n"
             "Args:\n"
             "    pdb (str): Path to PDB file\n"
             "    psf (str or List[str], optional): Path(s) to PSF file(s)\n"
             "    top (str, optional): Path to TOP file\n"
             "    itp (str or List[str], optional): Path(s) to ITP file(s)\n\n"
             "Note: You can combine different file types as needed.")
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