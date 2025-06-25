// src/bindings/io/IOParametersBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/IOModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace io {

void init_parameters_bindings(py::module& m, py::module& io_module) {
    // INPParser bindings
    auto inp_parser = py::class_<pygcmc::io::INPParser>(io_module, "INPParser")
        .def_static("parse_file", &pygcmc::io::INPParser::parse_file,
            py::arg("filename"),
            "Parse an input file and return a new Param object")
        .def_static("parse_string", &pygcmc::io::INPParser::parse_string,
            py::arg("content"),
            "Parse an input string and return a new Param object")
        .def_static("parse_to_param", &pygcmc::io::INPParser::parse_to_param,
            py::arg("filename"), py::arg("param"),
            "Parse an input file into an existing Param object")
        .def_static("parse_string_to_param", &pygcmc::io::INPParser::parse_string_to_param,
            py::arg("content"), py::arg("param"),
            "Parse an input string into an existing Param object");
    
    // Add to main module for backward compatibility
    m.attr("INPParser") = inp_parser;
}

} // namespace io
} // namespace bindings
} // namespace pygcmc