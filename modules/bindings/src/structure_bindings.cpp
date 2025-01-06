#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/project_residue.hpp"
#include "pygcmc/core/project_atom.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

void init_structure_bindings(py::module& m) {
    // Bind Structure class with proper wrapper support
    py::class_<Structure, std::shared_ptr<Structure>>(m, "Structure")
        .def(py::init<>())
        .def(py::init([](const std::string& pdb, const std::string& psf) {
                Structure structure;
                structure.read_pdb(pdb);
                structure.read_psf(psf);
                return structure;
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Create a structure and load from PDB and PSF files")
        .def(py::init([](const std::string& pdb, const std::vector<std::string>& psf) {
                return Structure::from_pdb_psf(pdb, psf);
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("psf"),
             "Create a structure and load from PDB and multiple PSF files")
        .def(py::init([](const std::string& pdb, const std::string& top) {
                Structure structure;
                structure.read_pdb(pdb);
                structure.read_top(top);
                return structure;
             }),
             py::kw_only(),
             py::arg("pdb"), py::arg("top"),
             "Create a structure and load from PDB and TOP files")
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
                
                return Structure::from_kwargs(cpp_kwargs);
             }),
             "Create a structure with flexible file loading options.\n\n"
             "Args:\n"
             "    pdb (str): Path to PDB file\n"
             "    psf (str or List[str], optional): Path(s) to PSF file(s)\n"
             "    top (str, optional): Path to TOP file\n"
             "    itp (str or List[str], optional): Path(s) to ITP file(s)\n\n"
             "Note: You can combine different file types as needed.")
        .def("apply_forcefield", &Structure::apply_forcefield)
        // Add structure loading methods
        .def("read_pdb_file", &Structure::read_pdb_file, "Read structure from PDB file")
        .def("read_top_file", &Structure::read_top_file, "Read topology from TOP file (with includes)")
        .def("read_top_file_without_includes", &Structure::read_top_file_without_includes, "Read topology from TOP file without includes")
        .def("read_top_file_with_includes", &Structure::read_top_file_with_includes, "Read topology from TOP file with includes (alias for read_top_file)")
        // Add PSF and ITP loading methods (both read_* and load_* variants)
        .def("read_psf", &Structure::read_psf, "Read topology from PSF file")
        .def("read_itp", &Structure::read_itp, "Read topology from ITP file")
        .def("load_psf", &Structure::load_psf, "Read topology from PSF file (alias for read_psf)")
        .def("load_itp", &Structure::load_itp, "Read topology from ITP file (alias for read_itp)")
        // Add alias methods for test compatibility
        .def("read_pdb", &Structure::read_pdb, "Read structure from PDB file (alias for read_pdb_file)")
        .def("read_top", &Structure::read_top, "Read topology from TOP file (alias for read_top_file)")
        .def("read_top_without_includes", &Structure::read_top_without_includes, "Read topology from TOP file without includes (alias)")
        // Add new PSF loading methods from System
        .def("load_structure_psf_auto", &Structure::load_structure_psf_auto,
             py::arg("pdb_file"), py::arg("psf_file"),
             "Load structure from PDB and PSF files with automatic detection of PSF type")
        .def("load_structure_psf_multi", &Structure::load_structure_psf_multi,
             py::arg("pdb_file"), py::arg("psf_file"),
             "Load structure from PDB and multi-residue PSF file")
        .def("load_structure_psf_single", &Structure::load_structure_psf_single,
             py::arg("pdb_file"), py::arg("psf_file"), py::arg("target_residue"),
             "Load structure from PDB and single-residue PSF file, applying to specified residue type")
        // Add static factory methods
        .def_static("from_pdb_psf", &Structure::from_pdb_psf,
             py::arg("pdb_file"), py::arg("psf_files"),
             "Create a structure from PDB and multiple PSF files")
        .def_static("from_pdb_psf_itp", &Structure::from_pdb_psf_itp,
             py::arg("pdb_file"), py::arg("psf_files"), py::arg("itp_file"),
             "Create a structure from PDB, multiple PSF files, and an ITP file")
        .def_static("from_pdb_psf_itps", &Structure::from_pdb_psf_itps,
             py::arg("pdb_file"), py::arg("psf_files"), py::arg("itp_files"),
             "Create a structure from PDB, multiple PSF files, and multiple ITP files")
        // Add PDB atom management methods
        .def("get_pdb_atom_count", &Structure::get_pdb_atom_count,
             "Get number of PDB atoms")
        .def("get_pdb_atom",
             py::overload_cast<size_t>(&Structure::get_pdb_atom, py::const_),
             py::arg("index"),
             "Get PDB atom by index (const)")
        .def("get_pdb_atom",
             py::overload_cast<size_t>(&Structure::get_pdb_atom),
             py::arg("index"),
             "Get PDB atom by index (mutable)")
        .def("add_pdb_atom", &Structure::add_pdb_atom,
             py::arg("atom"),
             "Add a PDB atom to the structure")
        .def("remove_pdb_atom", &Structure::remove_pdb_atom,
             py::arg("index"),
             "Remove a PDB atom by index")
        .def("get_pdb_atoms_by_residue", &Structure::get_pdb_atoms_by_residue,
             py::arg("residue_name"),
             "Get PDB atoms by residue name")
        .def("get_pdb_atoms_by_residue_sequence", &Structure::get_pdb_atoms_by_residue_sequence,
             py::arg("residue_name"), py::arg("sequence"),
             "Get PDB atoms by residue name and sequence number")
        .def("get_pdb_atoms_by_chain", &Structure::get_pdb_atoms_by_chain,
             py::arg("chain"),
             "Get PDB atoms by chain identifier")
        .def("clear_pdb_atoms", &Structure::clear_pdb_atoms,
             "Clear all PDB atoms")
        .def("has_pdb_atoms", &Structure::has_pdb_atoms,
             "Check if structure has any PDB atoms")
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
        .def("__len__", &Structure::get_num_atoms)
        .def("get_box", &Structure::get_box)
        .def("set_box", &Structure::set_box)
        // Add new data access methods
        .def("get_coordinates", [](Structure& self) {
            auto coords = self.get_coordinates();
            py::list result;
            for (const auto& coord : coords) {
                py::list xyz;
                xyz.append(coord[0]);
                xyz.append(coord[1]);
                xyz.append(coord[2]);
                result.append(xyz);
            }
            return result;
        })
        .def("get_box_vectors", [](Structure& self) {
            auto vectors = self.get_box_vectors();
            py::list result;
            for (const auto& vector : vectors) {
                py::list xyz;
                xyz.append(vector[0]);
                xyz.append(vector[1]);
                xyz.append(vector[2]);
                result.append(xyz);
            }
            return result;
        })
        .def("get_atoms_data", [](Structure& self) {
            py::list result;
            for (const auto& atom : self.atoms()) {
                if (atom) {
                    py::dict atom_data;
                    atom_data["serial"] = atom->serial;
                    atom_data["name"] = atom->name;
                    atom_data["residue"] = atom->residue;
                    atom_data["sequence"] = atom->sequence;
                    atom_data["chain"] = std::string(1, atom->chain);
                    atom_data["type"] = atom->type;
                    atom_data["topo_type"] = atom->topo_type;
                    atom_data["x"] = atom->x;
                    atom_data["y"] = atom->y;
                    atom_data["z"] = atom->z;
                    atom_data["topo_charge"] = atom->topo_charge;
                    atom_data["topo_mass"] = atom->topo_mass;
                    result.append(atom_data);
                }
            }
            return result;
        })
        .def("get_residues_data", [](Structure& self) {
            py::list result;
            for (const auto& residue : self.residues()) {
                if (residue) {
                    py::dict residue_data;
                    residue_data["name"] = residue->name;
                    residue_data["sequence_number"] = residue->sequence_number;
                    residue_data["chain_id"] = std::string(1, residue->chain_id);
                    residue_data["n_atoms"] = residue->atom_count();
                    
                    auto com = residue->center_of_mass();
                    py::list com_list;
                    com_list.append(com[0]);
                    com_list.append(com[1]);
                    com_list.append(com[2]);
                    residue_data["center_of_mass"] = com_list;
                    
                    result.append(residue_data);
                }
            }
            return result;
        })
        .def("get_energy_components", &Structure::get_energy_components)
        .def("get_atom_energy_contributions", [](Structure& self) {
            auto contributions = self.get_atom_energy_contributions();
            py::list result;
            
            for (const auto& [idx, vdw_e, elec_e, total_e] : contributions) {
                py::dict energy_data;
                energy_data["atom_index"] = idx;
                energy_data["vdw_energy"] = vdw_e;
                energy_data["electrostatic_energy"] = elec_e;
                energy_data["total_energy"] = total_e;
                result.append(energy_data);
            }
            
            return result;
        });
} 