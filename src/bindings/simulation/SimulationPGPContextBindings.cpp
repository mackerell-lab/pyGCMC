#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "platform/cpu/energy/pgp/PGPContext.hpp"
#include "model/montecarlo/MCMain.hpp"

namespace py = pybind11;
using namespace platform::cpu::energy::pgp;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_pgp_context_bindings(py::module& m) {
    // EnergyComponents struct
    py::class_<PGPContext::EnergyComponents>(m, "PGPEnergyComponents")
        .def_readwrite("total", &PGPContext::EnergyComponents::total)
        .def_readwrite("real_space", &PGPContext::EnergyComponents::real_space)
        .def_readwrite("reciprocal", &PGPContext::EnergyComponents::reciprocal)
        .def_readwrite("self", &PGPContext::EnergyComponents::self)
        .def_readwrite("vdw", &PGPContext::EnergyComponents::vdw)
        .def("__repr__", [](const PGPContext::EnergyComponents& e) {
            return "<PGPEnergyComponents total=" + std::to_string(e.total) + 
                   " real=" + std::to_string(e.real_space) +
                   " recip=" + std::to_string(e.reciprocal) +
                   " self=" + std::to_string(e.self) +
                   " vdw=" + std::to_string(e.vdw) + ">";
        });
    
    // PGPContext class
    py::class_<PGPContext>(m, "PGPContext")
        .def(py::init<>())
        .def("initialize", [](PGPContext& self,
                            double cutoff,
                            py::list box,
                            double alpha,
                            py::list meshSize,
                            double potential_cutoff,
                            py::list potentialGridSize,
                            int splineOrder,
                            double tolerance) {
            // Convert Python lists to arrays
            if (box.size() != 3) {
                throw std::runtime_error("Box must have 3 elements");
            }
            if (meshSize.size() != 3) {
                throw std::runtime_error("MeshSize must have 3 elements");
            }
            if (potentialGridSize.size() != 3) {
                throw std::runtime_error("PotentialGridSize must have 3 elements");
            }
            
            std::array<double, 3> box_arr = {
                box[0].cast<double>(),
                box[1].cast<double>(),
                box[2].cast<double>()
            };
            
            std::array<int, 3> mesh_arr = {
                meshSize[0].cast<int>(),
                meshSize[1].cast<int>(),
                meshSize[2].cast<int>()
            };
            
            std::array<int, 3> grid_arr = {
                potentialGridSize[0].cast<int>(),
                potentialGridSize[1].cast<int>(),
                potentialGridSize[2].cast<int>()
            };
            
            self.initialize(cutoff, box_arr, alpha, mesh_arr, 
                          potential_cutoff, grid_arr, splineOrder, tolerance);
        }, py::arg("cutoff"), py::arg("box"), py::arg("alpha"),
           py::arg("meshSize"), py::arg("potential_cutoff"),
           py::arg("potentialGridSize"), py::arg("splineOrder"),
           py::arg("tolerance"),
           "Initialize PGP parameters for this instance")
        
        .def("compute_system_energy", &PGPContext::computeSystemEnergy,
             py::arg("state"),
             "Compute total system energy using this PGP instance")
        
        .def("compute_movement_energy", &PGPContext::computeMovementEnergy,
             py::arg("state"), py::arg("movement_residues"),
             "Compute movement energy using this PGP instance")
        
        .def("precompute_grid_potential", &PGPContext::precomputeGridPotential,
             py::arg("state"), py::arg("atom_type"),
             "Precompute grid potential for a specific atom type");
    
    // Convenience function to create and return energy as dict
    m.def("compute_pgp_energy_with_context", [](PGPContext& context, pygcmc::model::montecarlo::MCState& state) -> py::dict {
        auto energy = context.computeSystemEnergy(state);
        py::dict result;
        result["total"] = energy.total;
        result["real_space"] = energy.real_space; 
        result["reciprocal"] = energy.reciprocal;
        result["self"] = energy.self;
        result["vdw"] = energy.vdw;
        return result;
    }, py::arg("context"), py::arg("state"),
       "Compute PGP energy and return as dictionary");
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc