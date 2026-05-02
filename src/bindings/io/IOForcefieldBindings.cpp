// src/bindings/io/IOForcefieldBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/IOModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace io {

void init_forcefield_bindings(py::module& m, py::module& io_module) {
    // PRMParser bindings
    auto prm_parser = py::class_<pygcmc::io::PRMParser>(m, "PRMParser")
        .def(py::init<>())
        .def_static("parse_string", &pygcmc::io::PRMParser::parse_string,
            py::arg("content"), py::arg("ff"))
        .def_static("parse_file", &pygcmc::io::PRMParser::parse_file,
            py::arg("filename"))
        .def_static("parse_files", &pygcmc::io::PRMParser::parse_files,
            py::arg("filenames"))
        .def_static("parse_file_to_forcefield", &pygcmc::io::PRMParser::parse_file_to_forcefield,
            py::arg("filename"), py::arg("ff"))
        .def("parse", &pygcmc::io::PRMParser::parse,
            py::arg("filename"), py::arg("ff"));

    // Add to both main module and io submodule for backward compatibility
    m.attr("PRMParser") = prm_parser;
    io_module.attr("PRMParser") = prm_parser;
}

} // namespace io
} // namespace bindings
} // namespace pygcmc
