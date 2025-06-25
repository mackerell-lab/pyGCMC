// src/bindings/model/ModelMontecarloBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace model {

void init_montecarlo_bindings(py::module& m, py::module&) {
    // Bind GCMCInfo
    py::class_<pygcmc::model::MCInfo>(m, "MCInfo")
        .def(py::init<>())
        .def_readwrite("mc_steps", &pygcmc::model::MCInfo::mcSteps)
        .def_property("box",
            [](const pygcmc::model::MCInfo& info) {
                return std::vector<float>{info.box[0], info.box[1], info.box[2]};
            },
            [](pygcmc::model::MCInfo& info, const std::vector<float>& box) {
                if (box.size() != 3) throw std::runtime_error("Box must have 3 dimensions");
                info.box[0] = box[0];
                info.box[1] = box[1];
                info.box[2] = box[2];
            })
        .def_readwrite("cutoff", &pygcmc::model::MCInfo::cutoff)
        .def_readwrite("beta", &pygcmc::model::MCInfo::beta)
        .def_readwrite("max_residues", &pygcmc::model::MCInfo::maxResidues)
        .def_readwrite("max_atoms", &pygcmc::model::MCInfo::maxAtoms)
        .def_readwrite("max_types", &pygcmc::model::MCInfo::maxTypes)
        .def_readwrite("volume", &pygcmc::model::MCInfo::volume)
        .def_readwrite("seed", &pygcmc::model::MCInfo::seed)
        .def_readwrite("use_switching", &pygcmc::model::MCInfo::use_switching)
        .def_readwrite("r_on", &pygcmc::model::MCInfo::r_on)
        .def_readwrite("r_off", &pygcmc::model::MCInfo::r_off)
        .def("setTemperature", &pygcmc::model::MCInfo::setTemperature, "Set temperature in Kelvin and calculate beta");

    // Bind GCMCInfo::Statistics
    py::class_<pygcmc::model::MCInfo::Statistics>(m, "MCStatistics")
        .def(py::init<>())
        .def_readwrite("totalMoves", &pygcmc::model::MCInfo::Statistics::totalMoves)
        .def_readwrite("acceptedMoves", &pygcmc::model::MCInfo::Statistics::acceptedMoves)
        .def_readwrite("insertionAttempts", &pygcmc::model::MCInfo::Statistics::insertionAttempts)
        .def_readwrite("acceptedInsertions", &pygcmc::model::MCInfo::Statistics::acceptedInsertions)
        .def_readwrite("deletionAttempts", &pygcmc::model::MCInfo::Statistics::deletionAttempts)
        .def_readwrite("acceptedDeletions", &pygcmc::model::MCInfo::Statistics::acceptedDeletions);

    // Bind MCResidue
    py::class_<pygcmc::model::MCResidue>(m, "MCResidue")
        .def(py::init<>())
        .def_readwrite("atomStart", &pygcmc::model::MCResidue::atomStart)
        .def_property_readonly("atom_start", [](const pygcmc::model::MCResidue& r) { return r.atomStart; })
        .def_readwrite("atomCount", &pygcmc::model::MCResidue::atomCount)
        .def_property_readonly("atom_count", [](const pygcmc::model::MCResidue& r) { return r.atomCount; })
        .def_readwrite("active", &pygcmc::model::MCResidue::active)
        .def_readwrite("fixed", &pygcmc::model::MCResidue::fixed)
        .def_property("center",
            [](const pygcmc::model::MCResidue& res) -> std::vector<float> {
                return {res.center[0], res.center[1], res.center[2]};
            },
            [](pygcmc::model::MCResidue& res, const std::vector<float>& center) {
                if (center.size() != 3) {
                    throw std::runtime_error("Center must be a vector of 3 floats");
                }
                res.center[0] = center[0];
                res.center[1] = center[1];
                res.center[2] = center[2];
            })
        .def_readwrite("concentration", &pygcmc::model::MCResidue::concentration)
        .def_readwrite("chemPot", &pygcmc::model::MCResidue::chemPot)
        .def_property_readonly("chem_pot", [](const pygcmc::model::MCResidue& r) { return r.chemPot; })
        .def_readwrite("type", &pygcmc::model::MCResidue::type)
        .def_readwrite("radius", &pygcmc::model::MCResidue::radius)
        .def_readwrite("energy_vdw", &pygcmc::model::MCResidue::energy_vdw)
        .def_readwrite("energy_elec", &pygcmc::model::MCResidue::energy_elec);

    // Bind MCAtom
    py::class_<pygcmc::model::MCAtom>(m, "MCAtom")
        .def(py::init<>())
        .def_readwrite("x", &pygcmc::model::MCAtom::x)
        .def_readwrite("y", &pygcmc::model::MCAtom::y)
        .def_readwrite("z", &pygcmc::model::MCAtom::z)
        .def_readwrite("charge", &pygcmc::model::MCAtom::charge)
        .def_readwrite("type", &pygcmc::model::MCAtom::type);

    // Bind MCMovementResidueInfo
    py::class_<pygcmc::model::MCMovementResidueInfo>(m, "MCMovementResidueInfo")
        .def(py::init<>())
        .def_readwrite("startIndex", &pygcmc::model::MCMovementResidueInfo::startIndex)
        .def_readwrite("activeCount", &pygcmc::model::MCMovementResidueInfo::activeCount)
        .def_readwrite("totalCount", &pygcmc::model::MCMovementResidueInfo::totalCount)
        .def_readwrite("resName", &pygcmc::model::MCMovementResidueInfo::resName);

    // Bind MCState
    py::class_<pygcmc::model::MCState>(m, "MCState")
        .def(py::init<>())
        .def("copy", [](const pygcmc::model::MCState& state) {
            pygcmc::model::MCState new_state;
            new_state.atoms = state.atoms;
            new_state.residues = state.residues;
            new_state.residueTypes = state.residueTypes;
            new_state.atomTypes = state.atomTypes;
            new_state.activeAtomCount = state.activeAtomCount;
            new_state.activeResidueCount = state.activeResidueCount;
            new_state.info = state.info;
            new_state.forcefield = state.forcefield;
            new_state.movementResidues = state.movementResidues;
            new_state.movementAtomTypes = state.movementAtomTypes;
            new_state.numMovementAtomTypes = state.numMovementAtomTypes;
            return new_state;
        }, "Create a deep copy of the MCState object")
        .def_property("atoms",
            [](pygcmc::model::MCState& state) -> std::vector<std::reference_wrapper<pygcmc::model::MCAtom>> {
                std::vector<std::reference_wrapper<pygcmc::model::MCAtom>> refs;
                refs.reserve(state.activeAtomCount);
                for (int i = 0; i < state.activeAtomCount; ++i) {
                    refs.push_back(std::ref(state.atoms[i]));
                }
                return refs;
            },
            [](pygcmc::model::MCState& state, const std::vector<pygcmc::model::MCAtom>& atoms) {
                state.atoms = atoms;
            })
        .def_property("residues",
            [](const pygcmc::model::MCState& state) {
                return std::vector<pygcmc::model::MCResidue>(state.residues.begin(), 
                    state.residues.begin() + state.activeResidueCount);
            },
            [](pygcmc::model::MCState& state, const std::vector<pygcmc::model::MCResidue>& residues) {
                state.residues = residues;
            })
        .def_readwrite("residueTypes", &pygcmc::model::MCState::residueTypes)
        .def_readwrite("atomTypes", &pygcmc::model::MCState::atomTypes)
        .def_readwrite("activeAtomCount", &pygcmc::model::MCState::activeAtomCount)
        .def_readwrite("activeResidueCount", &pygcmc::model::MCState::activeResidueCount)
        .def_readwrite("info", &pygcmc::model::MCState::info)
        .def_readwrite("forcefield", &pygcmc::model::MCState::forcefield)
        .def_property("movementResidues",
            [](const pygcmc::model::MCState& state) {
                return state.movementResidues;
            },
            [](pygcmc::model::MCState& state, const std::vector<pygcmc::model::MCMovementResidueInfo>& movementResidues) {
                state.movementResidues = movementResidues;
            })
        .def_readwrite("movementAtomTypes", &pygcmc::model::MCState::movementAtomTypes)
        .def_readwrite("numMovementAtomTypes", &pygcmc::model::MCState::numMovementAtomTypes)
        .def_property_readonly("ewald_energy", [](const pygcmc::model::MCState& state) {
            py::dict result;
            result["real_space"] = state.ewald_energy.real_space;
            result["reciprocal"] = state.ewald_energy.reciprocal;
            result["self"] = state.ewald_energy.self;
            result["total"] = state.ewald_energy.total;
            return result;
        });

    // Bind MCForceField
    py::class_<pygcmc::model::MCForceField>(m, "MCForceField")
        .def(py::init<>())
        .def_readwrite("numTotalTypes", &pygcmc::model::MCForceField::numTotalTypes)
        .def_property("maxTypes",
            [](const pygcmc::model::MCForceField& ff) { return ff.numTotalTypes; },
            [](pygcmc::model::MCForceField& ff, int value) { ff.numTotalTypes = value; })
        .def_readwrite("numMovementTypes", &pygcmc::model::MCForceField::numMovementTypes)
        .def_readwrite("ljSigma", &pygcmc::model::MCForceField::ljSigma)
        .def_readwrite("ljEps", &pygcmc::model::MCForceField::ljEps);

    // Add COULOMB constant to the module
    m.attr("COULOMB") = 138.935458; // kJ·mol^-1·nm·e^-2, Coulomb's constant in MD units
}

} // namespace model
} // namespace bindings
} // namespace pygcmc