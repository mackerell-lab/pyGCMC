#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pygcmc/core/system.hpp"
#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/io/psf_parser.hpp"

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
                // Create a structure
                Structure structure;

                // Load PDB file
                structure.read_pdb(pdb);

                // Read the first PSF file to establish base topology
                structure.read_psf(psf[0]);

                // Get all atom pointers after reading first PSF
                std::vector<io::PDBAtom*> atom_ptrs;
                for (const auto& residue : structure.residues()) {
                    for (const auto& atom : residue->atom_ptrs) {
                        if (atom) {
                            atom_ptrs.push_back(atom.get());
                        }
                    }
                }

                // Update atoms with remaining PSF files
                if (psf.size() > 1) {
                    std::vector<std::string> remaining_psf_files(psf.begin() + 1, psf.end());
                    int updated = io::PSFParser::update_pdb_atoms_from_multiple_psf(atom_ptrs, remaining_psf_files);
                    if (updated == 0) {
                        throw std::runtime_error("No atoms were updated with PSF information from additional PSF files");
                    }
                }

                // Create and return the system
                System system;
                system.load_structure(structure);
                return system;
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Create a system and load structure from PDB and multiple PSF files")
        .def(py::init([](const std::string& pdb, const std::vector<std::string>& psf, const std::string& itp) {
                // Create a structure
                Structure structure;

                // Load PDB file
                structure.read_pdb(pdb);

                // Read the first PSF file to establish base topology
                structure.read_psf(psf[0]);

                // Get all atom pointers after reading first PSF
                std::vector<io::PDBAtom*> atom_ptrs;
                for (const auto& residue : structure.residues()) {
                    for (const auto& atom : residue->atom_ptrs) {
                        if (atom) {
                            atom_ptrs.push_back(atom.get());
                        }
                    }
                }

                // Update atoms with remaining PSF files
                if (psf.size() > 1) {
                    std::vector<std::string> remaining_psf_files(psf.begin() + 1, psf.end());
                    int updated = io::PSFParser::update_pdb_atoms_from_multiple_psf(atom_ptrs, remaining_psf_files);
                    if (updated == 0) {
                        throw std::runtime_error("No atoms were updated with PSF information from additional PSF files");
                    }
                }

                // Load ITP file
                structure.read_itp(itp);

                // Create and return the system
                System system;
                system.load_structure(structure);
                return system;
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"), py::arg("itp"),
             "Create a system and load structure from PDB, multiple PSF files, and an ITP file")
        .def(py::init([](const std::string& pdb, const std::vector<std::string>& psf, const std::vector<std::string>& itp) {
                // Create a structure
                Structure structure;

                // Load PDB file
                structure.read_pdb(pdb);

                // Get all atom pointers
                std::vector<io::PDBAtom*> atom_ptrs;
                for (const auto& residue : structure.residues()) {
                    for (const auto& atom : residue->atom_ptrs) {
                        if (atom) {
                            atom_ptrs.push_back(atom.get());
                        }
                    }
                }

                // Try to update atoms using multiple PSF method
                int updated = io::PSFParser::update_pdb_atoms_from_multiple_psf(atom_ptrs, psf);
                if (updated == 0) {
                    throw std::runtime_error("No atoms were updated with PSF information using multiple PSF method");
                }

                // Load ITP files
                for (const auto& itp_file : itp) {
                    structure.read_itp(itp_file);
                }

                // Create and return the system
                System system;
                system.load_structure(structure);
                return system;
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
                 // Check for PDB file
                 if (!kwargs.contains("pdb")) {
                     throw std::runtime_error("PDB file is required");
                 }
                 std::string pdb_file = kwargs["pdb"].cast<std::string>();

                 // Handle PSF files
                 if (kwargs.contains("psf")) {
                     py::object psf_obj = kwargs["psf"];
                     if (py::isinstance<py::str>(psf_obj)) {
                         // Single PSF file
                         std::string psf_file = psf_obj.cast<std::string>();
                         self.load_structure_psf_auto(pdb_file, psf_file);
                     } else if (py::isinstance<py::list>(psf_obj)) {
                         // Multiple PSF files - use the same logic as constructor
                         std::vector<std::string> psf_files = psf_obj.cast<std::vector<std::string>>();
                         
                         // Create a structure
                         Structure structure;
                         
                         // Load PDB file
                         structure.read_pdb(pdb_file);
                         
                         // Read the first PSF file to establish base topology
                         structure.read_psf(psf_files[0]);
                         
                         // Get all atom pointers after reading first PSF
                         std::vector<io::PDBAtom*> atom_ptrs;
                         for (const auto& residue : structure.residues()) {
                             for (const auto& atom : residue->atom_ptrs) {
                                 if (atom) {
                                     atom_ptrs.push_back(atom.get());
                                 }
                             }
                         }
                         
                         // Update atoms with remaining PSF files
                         if (psf_files.size() > 1) {
                             std::vector<std::string> remaining_psf_files(psf_files.begin() + 1, psf_files.end());
                             int updated = io::PSFParser::update_pdb_atoms_from_multiple_psf(atom_ptrs, remaining_psf_files);
                             if (updated == 0) {
                                 throw std::runtime_error("No atoms were updated with PSF information from additional PSF files");
                             }
                         }
                         
                         // Load the structure into the system
                         self.load_structure(structure);
                     }
                 }

                 // Handle TOP file
                 if (kwargs.contains("top")) {
                     std::string top_file = kwargs["top"].cast<std::string>();
                     self.load_structure_top(pdb_file, top_file);
                 }

                 // Handle ITP files
                 if (kwargs.contains("itp")) {
                     py::object itp_obj = kwargs["itp"];
                     if (py::isinstance<py::str>(itp_obj)) {
                         // Single ITP file
                         std::string itp_file = itp_obj.cast<std::string>();
                         Structure structure;
                         structure.read_pdb(pdb_file);
                         structure.read_itp(itp_file);
                         self.load_structure(structure);
                     } else if (py::isinstance<py::list>(itp_obj)) {
                         // Multiple ITP files
                         Structure structure;
                         structure.read_pdb(pdb_file);
                         std::vector<std::string> itp_files = itp_obj.cast<std::vector<std::string>>();
                         for (const auto& itp_file : itp_files) {
                             structure.read_itp(itp_file);
                         }
                         self.load_structure(structure);
                     }
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
             "Check if system has any PDB atoms")
        .def(py::init([](py::kwargs kwargs) {
                // Create a structure
                Structure structure;
                System system;

                // Check for PDB file
                if (kwargs.contains("pdb")) {
                    std::string pdb_file = kwargs["pdb"].cast<std::string>();
                    structure.read_pdb(pdb_file);
                }

                // Get all atom pointers
                std::vector<io::PDBAtom*> atom_ptrs;
                for (const auto& residue : structure.residues()) {
                    for (const auto& atom : residue->atom_ptrs) {
                        if (atom) {
                            atom_ptrs.push_back(atom.get());
                        }
                    }
                }

                // Handle PSF files
                if (kwargs.contains("psf")) {
                    py::object psf_obj = kwargs["psf"];
                    if (py::isinstance<py::str>(psf_obj)) {
                        // Single PSF file
                        std::string psf_file = psf_obj.cast<std::string>();
                        structure.read_psf(psf_file);
                    } else if (py::isinstance<py::list>(psf_obj)) {
                        // Multiple PSF files
                        std::vector<std::string> psf_files = psf_obj.cast<std::vector<std::string>>();
                        int updated = io::PSFParser::update_pdb_atoms_from_multiple_psf(atom_ptrs, psf_files);
                        if (updated == 0) {
                            throw std::runtime_error("No atoms were updated with PSF information using multiple PSF method");
                        }
                    }
                }

                // Handle TOP file
                if (kwargs.contains("top")) {
                    std::string top_file = kwargs["top"].cast<std::string>();
                    structure.read_top(top_file);
                }

                // Handle ITP files
                if (kwargs.contains("itp")) {
                    py::object itp_obj = kwargs["itp"];
                    if (py::isinstance<py::str>(itp_obj)) {
                        // Single ITP file
                        std::string itp_file = itp_obj.cast<std::string>();
                        structure.read_itp(itp_file);
                    } else if (py::isinstance<py::list>(itp_obj)) {
                        // Multiple ITP files
                        std::vector<std::string> itp_files = itp_obj.cast<std::vector<std::string>>();
                        for (const auto& itp_file : itp_files) {
                            structure.read_itp(itp_file);
                        }
                    }
                }

                // Load the structure into the system
                system.load_structure(structure);
                return system;
             }),
             "Create a system with flexible file loading options.\n\n"
             "Args:\n"
             "    pdb (str): Path to PDB file\n"
             "    psf (str or List[str], optional): Path(s) to PSF file(s)\n"
             "    top (str, optional): Path to TOP file\n"
             "    itp (str or List[str], optional): Path(s) to ITP file(s)\n\n"
             "Note: You can combine different file types as needed.")
        .def("load_structure",
             [](System& self, py::kwargs kwargs) {
                 // Check for PDB file
                 if (!kwargs.contains("pdb")) {
                     throw std::runtime_error("PDB file is required");
                 }
                 std::string pdb_file = kwargs["pdb"].cast<std::string>();

                 // Handle PSF files
                 if (kwargs.contains("psf")) {
                     py::object psf_obj = kwargs["psf"];
                     if (py::isinstance<py::str>(psf_obj)) {
                         // Single PSF file
                         std::string psf_file = psf_obj.cast<std::string>();
                         self.load_structure_psf_auto(pdb_file, psf_file);
                     } else if (py::isinstance<py::list>(psf_obj)) {
                         // Multiple PSF files - use the same logic as constructor
                         std::vector<std::string> psf_files = psf_obj.cast<std::vector<std::string>>();
                         
                         // Create a structure
                         Structure structure;
                         
                         // Load PDB file
                         structure.read_pdb(pdb_file);
                         
                         // Read the first PSF file to establish base topology
                         structure.read_psf(psf_files[0]);
                         
                         // Get all atom pointers after reading first PSF
                         std::vector<io::PDBAtom*> atom_ptrs;
                         for (const auto& residue : structure.residues()) {
                             for (const auto& atom : residue->atom_ptrs) {
                                 if (atom) {
                                     atom_ptrs.push_back(atom.get());
                                 }
                             }
                         }
                         
                         // Update atoms with remaining PSF files
                         if (psf_files.size() > 1) {
                             std::vector<std::string> remaining_psf_files(psf_files.begin() + 1, psf_files.end());
                             int updated = io::PSFParser::update_pdb_atoms_from_multiple_psf(atom_ptrs, remaining_psf_files);
                             if (updated == 0) {
                                 throw std::runtime_error("No atoms were updated with PSF information from additional PSF files");
                             }
                         }
                         
                         // Load the structure into the system
                         self.load_structure(structure);
                     }
                 }

                 // Handle TOP file
                 if (kwargs.contains("top")) {
                     std::string top_file = kwargs["top"].cast<std::string>();
                     self.load_structure_top(pdb_file, top_file);
                 }

                 // Handle ITP files
                 if (kwargs.contains("itp")) {
                     py::object itp_obj = kwargs["itp"];
                     if (py::isinstance<py::str>(itp_obj)) {
                         // Single ITP file
                         std::string itp_file = itp_obj.cast<std::string>();
                         Structure structure;
                         structure.read_pdb(pdb_file);
                         structure.read_itp(itp_file);
                         self.load_structure(structure);
                     } else if (py::isinstance<py::list>(itp_obj)) {
                         // Multiple ITP files
                         Structure structure;
                         structure.read_pdb(pdb_file);
                         std::vector<std::string> itp_files = itp_obj.cast<std::vector<std::string>>();
                         for (const auto& itp_file : itp_files) {
                             structure.read_itp(itp_file);
                         }
                         self.load_structure(structure);
                     }
                 }
             },
             "Load structure with flexible file loading options.\n\n"
             "Args:\n"
             "    pdb (str): Path to PDB file\n"
             "    psf (str or List[str], optional): Path(s) to PSF file(s)\n"
             "    top (str, optional): Path to TOP file\n"
             "    itp (str or List[str], optional): Path(s) to ITP file(s)\n\n"
             "Note: You can combine different file types as needed.");
} 