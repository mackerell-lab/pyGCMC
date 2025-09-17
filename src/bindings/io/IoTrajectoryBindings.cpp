// src/bindings/io/IoTrajectoryBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../io/output/TrajectoryWriter.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace io {

void init_trajectory_bindings(py::module& m) {
    
    // Create a wrapper class that provides the old interface
    class TrajectoryWriterWrapper {
    public:
        enum class Format {
            PDB,
            XYZ,
            DAT,
            TOP
        };
        
    private:
        std::unique_ptr<pygcmc::io::output::TrajectoryWriter> writer_;
        Format format_;
        std::string filename_;
        
    public:
        TrajectoryWriterWrapper(const std::string& filename, Format format = Format::PDB) 
            : format_(format), filename_(filename) {
            pygcmc::io::output::TrajectoryWriter::Config config;
            
            switch(format) {
                case Format::PDB:
                    config.format = "pdb";
                    break;
                case Format::XYZ:
                    config.format = "xyz";
                    break;
                case Format::DAT:
                    config.format = "dat";
                    break;
                case Format::TOP:
                    config.format = "top";
                    break;
            }
            
            config.continuousFile = true;  // Keep file open for compatibility
            config.multiFrame = (format == Format::PDB);  // PDB supports multi-frame
            
            writer_ = std::make_unique<pygcmc::io::output::TrajectoryWriter>(config);
            writer_->open(filename);
        }
        
        void writeFrame(const model::montecarlo::MCState& state, int frameNumber = -1) {
            if (writer_) {
                writer_->writeTrajectory(state, frameNumber, "");
            }
        }
        
        void writeStatistics(const model::montecarlo::MCInfo::Statistics& /* stats */) {
            // Not directly supported in new interface, ignore for compatibility
        }
        
        void writeTopology(const model::montecarlo::MCState& state) {
            if (writer_ && format_ == Format::TOP) {
                writer_->writeTrajectory(state, 0, "");
            }
        }
        
        void close() {
            if (writer_) {
                writer_->close();
            }
        }
        
        // Context manager support
        TrajectoryWriterWrapper& enter() { return *this; }
        void exit(py::object, py::object, py::object) { close(); }
    };
    
    // Bind the wrapper class as TrajectoryWriter for backward compatibility
    py::class_<TrajectoryWriterWrapper> writer(m, "TrajectoryWriter");
    
    // Bind the Format enum
    py::enum_<TrajectoryWriterWrapper::Format>(writer, "Format")
        .value("PDB", TrajectoryWriterWrapper::Format::PDB)
        .value("XYZ", TrajectoryWriterWrapper::Format::XYZ)
        .value("DAT", TrajectoryWriterWrapper::Format::DAT)
        .value("TOP", TrajectoryWriterWrapper::Format::TOP);
    
    writer
        .def(py::init<const std::string&, TrajectoryWriterWrapper::Format>(),
             py::arg("filename"),
             py::arg("format") = TrajectoryWriterWrapper::Format::PDB,
             "Create a trajectory writer")
        .def(py::init<const std::string&>(),
             py::arg("filename"),
             "Create a trajectory writer with default PDB format")
        .def("write_frame", &TrajectoryWriterWrapper::writeFrame,
             py::arg("state"),
             py::arg("frame_number") = -1,
             "Write a single frame to the trajectory")
        .def("writeFrame", &TrajectoryWriterWrapper::writeFrame,
             py::arg("state"),
             py::arg("frame_number") = -1,
             "Write a single frame to the trajectory")
        .def("write_statistics", &TrajectoryWriterWrapper::writeStatistics,
             py::arg("stats"),
             "Write statistics data")
        .def("writeStatistics", &TrajectoryWriterWrapper::writeStatistics,
             py::arg("stats"),
             "Write statistics data")
        .def("write_topology", &TrajectoryWriterWrapper::writeTopology,
             py::arg("state"),
             "Write topology information")
        .def("writeTopology", &TrajectoryWriterWrapper::writeTopology,
             py::arg("state"),
             "Write topology information")
        .def("close", &TrajectoryWriterWrapper::close,
             "Close the trajectory file")
        .def("__enter__", &TrajectoryWriterWrapper::enter,
             py::return_value_policy::reference_internal)
        .def("__exit__", &TrajectoryWriterWrapper::exit);
    
    // Also bind the new interface
    py::class_<pygcmc::io::output::TrajectoryWriter::Config>(m, "TrajectoryWriterConfig")
        .def(py::init<>())
        .def_readwrite("format", &pygcmc::io::output::TrajectoryWriter::Config::format)
        .def_readwrite("prefix", &pygcmc::io::output::TrajectoryWriter::Config::prefix)
        .def_readwrite("compress_output", &pygcmc::io::output::TrajectoryWriter::Config::compressOutput)
        .def_readwrite("precision", &pygcmc::io::output::TrajectoryWriter::Config::precision)
        .def_readwrite("multi_frame", &pygcmc::io::output::TrajectoryWriter::Config::multiFrame)
        .def_readwrite("continuous_file", &pygcmc::io::output::TrajectoryWriter::Config::continuousFile);
    
    py::class_<pygcmc::io::output::TrajectoryWriter>(m, "TrajectoryWriterNew")
        .def(py::init<const pygcmc::io::output::TrajectoryWriter::Config&>(),
             py::arg("config"),
             "Create a trajectory writer with specified configuration")
        .def("open", &pygcmc::io::output::TrajectoryWriter::open,
             py::arg("filename") = "",
             "Open a trajectory file")
        .def("close", &pygcmc::io::output::TrajectoryWriter::close,
             "Close the trajectory file")
        .def("is_open", &pygcmc::io::output::TrajectoryWriter::isOpen,
             "Check if file is open")
        .def("write_trajectory", &pygcmc::io::output::TrajectoryWriter::writeTrajectory,
             py::arg("state"),
             py::arg("step"),
             py::arg("filename") = "",
             "Write a trajectory frame")
        .def("set_config", &pygcmc::io::output::TrajectoryWriter::setConfig,
             py::arg("config"),
             "Update the writer configuration")
        .def("get_config", &pygcmc::io::output::TrajectoryWriter::getConfig,
             "Get the current configuration")
        .def("get_frame_count", &pygcmc::io::output::TrajectoryWriter::getFrameCount,
             "Get the number of frames written")
        .def("reset_frame_count", &pygcmc::io::output::TrajectoryWriter::resetFrameCount,
             "Reset the frame counter");
    
    // Bind DataWriter
    py::class_<pygcmc::io::output::DataWriter>(m, "DataWriter")
        .def(py::init<const std::string&>(),
             py::arg("filename"),
             "Create a data writer for analysis output")
        .def("write_header", &pygcmc::io::output::DataWriter::writeHeader,
             py::arg("columns"),
             "Write column headers")
        .def("writeHeader", &pygcmc::io::output::DataWriter::writeHeader,
             py::arg("columns"),
             "Write column headers")
        .def("write_row", &pygcmc::io::output::DataWriter::writeRow,
             py::arg("values"),
             "Write a row of data")
        .def("writeRow", &pygcmc::io::output::DataWriter::writeRow,
             py::arg("values"),
             "Write a row of data")
        .def("write_comment", &pygcmc::io::output::DataWriter::writeComment,
             py::arg("comment"),
             "Write a comment line")
        .def("writeComment", &pygcmc::io::output::DataWriter::writeComment,
             py::arg("comment"),
             "Write a comment line")
        .def("close", &pygcmc::io::output::DataWriter::close,
             "Close the data file")
        .def("__enter__", [](pygcmc::io::output::DataWriter& self) -> pygcmc::io::output::DataWriter& {
            return self;
        }, py::return_value_policy::reference_internal)
        .def("__exit__", [](pygcmc::io::output::DataWriter& self, py::object, py::object, py::object) {
            self.close();
        });
}

} // namespace io
} // namespace bindings
} // namespace pygcmc