// src/bindings/io/IOTopologyBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/IOModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace io {

void init_topology_bindings(py::module& m, py::module& io_module) {
    // PSFParser bindings
    auto psf_parser = py::class_<pygcmc::io::PSFParser>(m, "PSFParser")
        .def(py::init<>())
        .def("parse_to_topology", &pygcmc::io::PSFParser::parse_to_topology)
        .def_static("parse_file", &pygcmc::io::PSFParser::parse_file,
            py::arg("filename"),
            "Parse a PSF file and return a new Topology object")
        .def_static("parse_string", &pygcmc::io::PSFParser::parse_string,
            py::arg("psf_str"),
            "Parse a PSF string and return a new Topology object");
    io_module.attr("PSFParser") = psf_parser;

    // TOPParser bindings
    auto top_parser = py::class_<pygcmc::io::TOPParser>(m, "TOPParser")
        .def(py::init<>())
        .def("parse_to_topology", &pygcmc::io::TOPParser::parse_to_topology)
        .def_static("parse_file", &pygcmc::io::TOPParser::parse_file,
            py::arg("filename"),
            "Parse a topology file and return a new Topology object")
        .def_static("parse_string", &pygcmc::io::TOPParser::parse_string,
            py::arg("top_str"),
            "Parse a topology string and return a new Topology object")
        .def_static("enable_debug", &pygcmc::io::TOPParser::enable_debug,
            py::arg("enable"),
            "Enable or disable debug output")
        .def_static("is_debug_enabled", &pygcmc::io::TOPParser::is_debug_enabled,
            "Check if debug output is enabled");
    io_module.attr("TOPParser") = top_parser;
}

} // namespace io
} // namespace bindings
} // namespace pygcmc