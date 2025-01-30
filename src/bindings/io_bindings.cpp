// src/bindings/io_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/pdbParser.hpp"
#include "model/structure.hpp"
#include "io/psfParser.hpp"
#include "io/topParser.hpp"
#include "model/atom.hpp"
#include "model/residue.hpp"
#include "../model/forcefield.hpp"
#include "../io/prmParser.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_io(py::module& m) {
    // Create io submodule
    auto io = m.def_submodule("io", "Input/Output operations");
    
    // Bind Structure::SecondaryStructure
    auto secondary_structure = py::class_<model::Structure::SecondaryStructure>(m, "SecondaryStructure")
        .def(py::init<>())
        .def_readwrite("id", &model::Structure::SecondaryStructure::id)
        .def_readwrite("initResName", &model::Structure::SecondaryStructure::initResName)
        .def_readwrite("initChainId", &model::Structure::SecondaryStructure::initChainId)
        .def_readwrite("initSeqNum", &model::Structure::SecondaryStructure::initSeqNum)
        .def_readwrite("initICode", &model::Structure::SecondaryStructure::initICode)
        .def_readwrite("endResName", &model::Structure::SecondaryStructure::endResName)
        .def_readwrite("endChainId", &model::Structure::SecondaryStructure::endChainId)
        .def_readwrite("endSeqNum", &model::Structure::SecondaryStructure::endSeqNum)
        .def_readwrite("endICode", &model::Structure::SecondaryStructure::endICode)
        .def_readwrite("helixClass", &model::Structure::SecondaryStructure::structureClass);

    // Bind Structure::TerminalInfo
    auto terminal_info = py::class_<model::Structure::TerminalInfo>(m, "TerminalInfo")
        .def(py::init<>())
        .def_readwrite("chainId", &model::Structure::TerminalInfo::chainId)
        .def_readwrite("resSeq", &model::Structure::TerminalInfo::resSeq)
        .def_readwrite("iCode", &model::Structure::TerminalInfo::iCode)
        .def_readwrite("resName", &model::Structure::TerminalInfo::resName);

    // Bind Structure to main module
    auto structure = py::class_<model::Structure>(m, "Structure")
        .def(py::init<>())
        // 直接暴露内部成员作为属性
        .def_property_readonly("atoms", [](const model::Structure& s) { return s.getAtoms(); },
            "List of atoms in the structure")
        .def_property_readonly("residues", [](const model::Structure& s) { return s.getResidues(); },
            "List of residues in the structure")
        .def_property_readonly("terminals", [](const model::Structure& s) { return s.getTerminals(); },
            "List of terminal records")
        .def_property_readonly("helices", [](const model::Structure& s) { return s.getHelices(); },
            "Map of chain IDs to helix information")
        .def_property_readonly("sheets", [](const model::Structure& s) { return s.getSheets(); },
            "Map of chain IDs to sheet information")
        .def_property_readonly("ssbonds", [](const model::Structure& s) { return s.getSSBonds(); },
            "List of disulfide bonds")
        .def_property_readonly("boxDimensions", [](const model::Structure& s) { return s.getBoxDimensions(); },
            "Box dimensions and angles")
        // 保留原有方法
        .def("addAtom", &model::Structure::addAtom)
        .def("addResidue", &model::Structure::addResidue)
        .def("addTerminal", &model::Structure::addTerminal)
        .def("addHelix", &model::Structure::addHelix)
        .def("addSheet", &model::Structure::addSheet)
        .def("addSSBond", &model::Structure::addSSBond)
        .def("setBoxDimensions", &model::Structure::setBoxDimensions)
        .def("clear", &model::Structure::clear);

    // Bind PDBParser to io submodule
    auto pdb_parser = py::class_<io::PDBParser>(io, "PDBParser")
        .def_static("parse_file", &io::PDBParser::parse_file,
            py::arg("filename"),
            "Parse PDB file and return Structure object")
        .def_static("parse_string", &io::PDBParser::parse_string,
            py::arg("pdbStr"),
            "Parse PDB string and return Structure object")
        .def_static("parse_to_structure", &io::PDBParser::parse_to_structure,
            py::arg("filename"), py::arg("structure"),
            "Parse PDB file and populate Structure object")
        .def_static("parse_string_to_structure", &io::PDBParser::parse_string_to_structure,
            py::arg("pdbStr"), py::arg("structure"),
            "Parse PDB string and populate Structure object");
    m.attr("PDBParser") = pdb_parser;  // 将 PDBParser 也添加到主模块

    // PSFParser bindings
    auto psf_parser = py::class_<io::PSFParser>(m, "PSFParser")
        .def(py::init<>())
        .def("parse_to_topology", &io::PSFParser::parse_to_topology)
        .def_static("parse_file", &io::PSFParser::parse_file,
            py::arg("filename"),
            "Parse a PSF file and return a new Topology object")
        .def_static("parse_string", &io::PSFParser::parse_string,
            py::arg("psf_str"),
            "Parse a PSF string and return a new Topology object");
    io.attr("PSFParser") = psf_parser;

    // TOPParser bindings
    auto top_parser = py::class_<io::TOPParser>(m, "TOPParser")
        .def(py::init<>())
        .def("parse_to_topology", &io::TOPParser::parse_to_topology)
        .def_static("parse_file", &io::TOPParser::parse_file,
            py::arg("filename"),
            "Parse a topology file and return a new Topology object")
        .def_static("parse_string", &io::TOPParser::parse_string,
            py::arg("top_str"),
            "Parse a topology string and return a new Topology object")
        .def_static("enable_debug", &io::TOPParser::enable_debug,
            py::arg("enable"),
            "Enable or disable debug output")
        .def_static("is_debug_enabled", &io::TOPParser::is_debug_enabled,
            "Check if debug output is enabled");
    io.attr("TOPParser") = top_parser;

    // Add PRMParser bindings
    auto prm_parser = py::class_<PRMParser>(m, "PRMParser")
        // Static methods
        .def_static("parse_string", &PRMParser::parse_string,
            "Parse parameter content from a string",
            py::arg("content"), py::arg("ff"))
        .def_static("parse_file", &PRMParser::parse_file,
            "Parse parameter content from a file",
            py::arg("filename"), py::arg("ff"))
        // Instance methods
        .def("parse", &PRMParser::parse,
            "Parse parameter content from a file (instance method)",
            py::arg("filename"), py::arg("ff"));

    m.attr("PRMParser") = prm_parser;  // Add to main module as well
    io.attr("PRMParser") = prm_parser;  // Add to io submodule
}

void init_forcefield(py::module& m) {
    // NonbondedParams
    py::class_<NonbondedParams>(m, "NonbondedParams")
        .def(py::init<>())
        .def_readwrite("nbxmod", &NonbondedParams::nbxmod)
        .def_readwrite("cdiel", &NonbondedParams::cdiel)
        .def_readwrite("fshift", &NonbondedParams::fshift)
        .def_readwrite("vatom", &NonbondedParams::vatom)
        .def_readwrite("vdistance", &NonbondedParams::vdistance)
        .def_readwrite("vfswitch", &NonbondedParams::vfswitch)
        .def_readwrite("cutnb", &NonbondedParams::cutnb)
        .def_readwrite("ctofnb", &NonbondedParams::ctofnb)
        .def_readwrite("ctonnb", &NonbondedParams::ctonnb)
        .def_readwrite("eps", &NonbondedParams::eps)
        .def_readwrite("e14fac", &NonbondedParams::e14fac)
        .def_readwrite("wmin", &NonbondedParams::wmin);

    // LJParams
    py::class_<LJParams>(m, "LJParams")
        .def(py::init<>())
        .def_readwrite("epsilon", &LJParams::epsilon)
        .def_readwrite("rmin", &LJParams::rmin);

    // BondParams
    py::class_<BondParams>(m, "BondParams")
        .def(py::init<>())
        .def_readwrite("kb", &BondParams::kb)
        .def_readwrite("b0", &BondParams::b0);

    // AngleParams
    py::class_<AngleParams>(m, "AngleParams")
        .def(py::init<>())
        .def_readwrite("ktheta", &AngleParams::ktheta)
        .def_readwrite("theta0", &AngleParams::theta0)
        .def_readwrite("kub", &AngleParams::kub)
        .def_readwrite("s0", &AngleParams::s0);

    // DihedralParams
    py::class_<DihedralParams>(m, "DihedralParams")
        .def(py::init<>())
        .def_readwrite("kchi", &DihedralParams::kchi)
        .def_readwrite("n", &DihedralParams::n)
        .def_readwrite("delta", &DihedralParams::delta);

    // ImproperParams
    py::class_<ImproperParams>(m, "ImproperParams")
        .def(py::init<>())
        .def_readwrite("kpsi", &ImproperParams::kpsi)
        .def_readwrite("psi0", &ImproperParams::psi0);

    // ForceField
    py::class_<ForceField>(m, "ForceField")
        .def(py::init<>())
        .def_readwrite("atom_masses", &ForceField::atom_masses)
        .def_readwrite("lj_params", &ForceField::lj_params)
        .def_readwrite("nbfix", &ForceField::nbfix)
        .def_readwrite("bond_params", &ForceField::bond_params)
        .def_readwrite("angle_params", &ForceField::angle_params)
        .def_readwrite("dihedral_params", &ForceField::dihedral_params)
        .def_readwrite("improper_params", &ForceField::improper_params)
        .def_readwrite("nonbonded_params", &ForceField::nonbonded_params)
        .def("get_nonbonded_params", [](const ForceField& self) -> const NonbondedParams& {
            return self.nonbonded_params;
        })
        .def("get_lj_params", [](const ForceField& self, const std::string& type) -> const LJParams& {
            auto it = self.lj_params.find(type);
            if (it == self.lj_params.end()) {
                throw py::key_error("No LJ parameters found for atom type: " + type);
            }
            return it->second;
        })
        .def("get_nbfix", [](const ForceField& self, const std::string& type1, const std::string& type2) -> py::tuple {
            auto key = std::make_pair(type1 < type2 ? type1 : type2, type1 < type2 ? type2 : type1);
            auto it = self.nbfix.find(key);
            if (it == self.nbfix.end()) {
                return py::make_tuple(0.0, false);
            }
            return py::make_tuple(it->second, true);
        });
}

} // namespace bindings
} // namespace pygcmc


