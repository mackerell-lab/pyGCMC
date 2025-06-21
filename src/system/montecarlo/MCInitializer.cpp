#include "MCInitializer.hpp" 
#include "MCCore.hpp"
#include <stdexcept>
#include <set>

namespace pygcmc {
namespace system {
namespace montecarlo {

void MCInitializer::initializeFromMolecular(model::MCState& state, const std::shared_ptr<model::Molecular>& molecular) {
    if (!molecular) {
        throw std::runtime_error("MolecularSystem has no molecular data");
    }
    
    // Set box dimensions from molecular system
    state.info.box[0] = molecular->boxDimensions[0] * ANGSTROM_TO_NM;  // Convert Å to nm
    state.info.box[1] = molecular->boxDimensions[1] * ANGSTROM_TO_NM;  // Convert Å to nm
    state.info.box[2] = molecular->boxDimensions[2] * ANGSTROM_TO_NM;  // Convert Å to nm
    state.info.volume = state.info.box[0] * state.info.box[1] * state.info.box[2];

    // Convert residues and atoms
    auto convertedData = convertMolecularData(state, molecular);
    
    // Check capacity
    if (convertedData.first.size() > static_cast<size_t>(state.info.maxResidues) ||
        convertedData.second.size() > static_cast<size_t>(state.info.maxAtoms)) {
        throw std::runtime_error("Initial system exceeds max capacity");
    }

    // Initialize the system with converted data using MCCore
    MCCore core;
    core.addInitialResidues(state, convertedData.first.data(), convertedData.first.size(),
                           convertedData.second.data(), convertedData.second.size());
}

void MCInitializer::initializeForceField(model::MCState& state, const model::ForceField& ff) {
    // TODO: Implement force field initialization
    // This would involve:
    // 1. Unit conversions (Å -> nm, kcal/mol -> kJ/mol)
    // 2. Transform Rmin/2 to sigma
    // 3. Apply combining rules or NBFIX
    // 4. Organize parameters for efficient access
    (void)state; // Suppress unused warning for now
    (void)ff;    // Suppress unused warning for now
}

void MCInitializer::validateParameters(const model::ForceField& ff, const std::shared_ptr<model::Molecular>& molecular) {
    if (!molecular) {
        throw std::runtime_error("Molecular system is null");
    }

    // Check all atom topology parameters
    for (const auto& res : molecular->topology_residues) {
        for (int atomIdx : res.atoms) {
            if (atomIdx >= static_cast<int>(molecular->topology_atoms.size())) {
                throw std::runtime_error("Invalid topology atom index: " + std::to_string(atomIdx));
            }
        }
    }

    // Check all atom type LJ parameters
    std::set<std::string> missingLJTypes;
    for (const auto& atom : molecular->topology_atoms) {
        try {
            ff.get_lj_params(atom.type);
        } catch (const std::exception&) {
            missingLJTypes.insert(atom.type);
        }
    }

    if (!missingLJTypes.empty()) {
        std::string error = "Missing LJ parameters for atom types: ";
        for (const auto& type : missingLJTypes) {
            error += type + " ";
        }
        throw std::runtime_error(error);
    }
}

std::pair<std::vector<model::MCResidue>, std::vector<model::MCAtom>>
MCInitializer::convertMolecularData(model::MCState& state, const std::shared_ptr<model::Molecular>& molecular) {
    std::vector<model::MCResidue> tempResidues;
    std::vector<model::MCAtom> tempAtoms;
    
    size_t atomStart = 0;
    const size_t numResidues = molecular->get_num_residues();
    
    for (size_t i = 0; i < numResidues; ++i) {
        const auto& molRes = molecular->residues[i];
        const auto& topRes = molecular->topology_residues[i];
        
        model::MCResidue mcRes;
        mcRes.atomStart = atomStart;
        mcRes.atomCount = molRes->atom_count();
        mcRes.active = true;
        
        // Initialize energy components and GCMC parameters in GROMACS units
        mcRes.energy_vdw = 0.0f;   // kJ/mole
        mcRes.energy_elec = 0.0f;  // kJ/mole
        mcRes.chemPot = 0.0f;      // kJ/mole
        mcRes.concentration = 0.0f; // mol/L
        mcRes.radius = 0.0f;       // nm
        
        // Convert atoms for this residue
        const auto& molAtoms = molRes->get_atoms();
        for (size_t j = 0; j < molAtoms.size(); j++) {
            const auto& molAtom = molAtoms[j];
            const auto& topAtom = molecular->topology_atoms[topRes.atoms[j]];
            
            model::MCAtom mcAtom;
            // Convert coordinates from Å to nm
            mcAtom.x = molAtom->get_x() * ANGSTROM_TO_NM;
            mcAtom.y = molAtom->get_y() * ANGSTROM_TO_NM;
            mcAtom.z = molAtom->get_z() * ANGSTROM_TO_NM;
            mcAtom.charge = topAtom.charge;  // Charge unit (e) remains the same
            mcAtom.type = state.atomTypes.getOrAddType(topAtom.type);
            
            tempAtoms.push_back(mcAtom);
        }

        // Set residue type
        mcRes.type = state.residueTypes.getOrAddType(molRes->get_resname());

        // Calculate center of mass (in nm)
        mcRes.center[0] = mcRes.center[1] = mcRes.center[2] = 0.0f;
        for (const auto& atom : molAtoms) {
            mcRes.center[0] += atom->get_x() * ANGSTROM_TO_NM;
            mcRes.center[1] += atom->get_y() * ANGSTROM_TO_NM;
            mcRes.center[2] += atom->get_z() * ANGSTROM_TO_NM;
        }
        
        if (mcRes.atomCount > 0) {
            float invCount = 1.0f / mcRes.atomCount;
            mcRes.center[0] *= invCount;
            mcRes.center[1] *= invCount;
            mcRes.center[2] *= invCount;
        }
        
        tempResidues.push_back(mcRes);
        atomStart += mcRes.atomCount;
    }

    return std::make_pair(std::move(tempResidues), std::move(tempAtoms));
}

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 