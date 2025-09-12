// src/bindings/io/IoTrajectoryBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/TrajectoryWriter.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace io {

void init_trajectory_bindings(py::module& m) {
    
    // Bind TrajectoryWriter
    py::class_<pygcmc::io::TrajectoryWriter> writer(m, "TrajectoryWriter");
    
    py::enum_<pygcmc::io::TrajectoryWriter::Format>(writer, "Format")
        .value("PDB", pygcmc::io::TrajectoryWriter::Format::PDB)
        .value("XYZ", pygcmc::io::TrajectoryWriter::Format::XYZ)
        .value("DAT", pygcmc::io::TrajectoryWriter::Format::DAT)
        .value("TOP", pygcmc::io::TrajectoryWriter::Format::TOP);
    
    writer
        .def(py::init<const std::string&, pygcmc::io::TrajectoryWriter::Format>(),
             py::arg("filename"),
             py::arg("format") = pygcmc::io::TrajectoryWriter::Format::PDB,
             "Create a trajectory writer")
        .def("write_frame", &pygcmc::io::TrajectoryWriter::writeFrame,
             py::arg("state"),
             py::arg("frame_number") = -1,
             "Write a single frame to the trajectory")
        .def("write_statistics", &pygcmc::io::TrajectoryWriter::writeStatistics,
             py::arg("stats"),
             "Write statistics data")
        .def("write_topology", &pygcmc::io::TrajectoryWriter::writeTopology,
             py::arg("state"),
             "Write topology information")
        .def("close", &pygcmc::io::TrajectoryWriter::close,
             "Close the trajectory file")
        .def("__enter__", [](pygcmc::io::TrajectoryWriter& self) -> pygcmc::io::TrajectoryWriter& {
            return self;
        })
        .def("__exit__", [](pygcmc::io::TrajectoryWriter& self, py::object, py::object, py::object) {
            self.close();
        });
    
    // Bind DataWriter
    py::class_<pygcmc::io::DataWriter>(m, "DataWriter")
        .def(py::init<const std::string&>(),
             py::arg("filename"),
             "Create a data writer for analysis output")
        .def("write_header", &pygcmc::io::DataWriter::writeHeader,
             py::arg("columns"),
             "Write column headers")
        .def("write_row", &pygcmc::io::DataWriter::writeRow,
             py::arg("values"),
             "Write a row of data")
        .def("write_comment", &pygcmc::io::DataWriter::writeComment,
             py::arg("comment"),
             "Write a comment line")
        .def("close", &pygcmc::io::DataWriter::close,
             "Close the data file")
        .def("__enter__", [](pygcmc::io::DataWriter& self) -> pygcmc::io::DataWriter& {
            return self;
        })
        .def("__exit__", [](pygcmc::io::DataWriter& self, py::object, py::object, py::object) {
            self.close();
        });
}

} // namespace io
} // namespace bindings
} // namespace pygcmc