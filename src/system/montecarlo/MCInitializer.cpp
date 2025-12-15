#include "MCInitializer.hpp" 
#include "MCCore.hpp"
#include <stdexcept>
#include <set>
#include <cmath>

namespace pygcmc {
namespace system {
namespace montecarlo {

void MCInitializer::initializeFromMolecular(model::MCState& state, const std::shared_ptr<model::Molecular>& molecular) {
    if (!molecular) {
        throw std::runtime_error("MolecularSystem has no molecular data");
    }
    
    // Set box dimensions from molecular system
    if (molecular->boxDimensions.size() >= 3) {
        state.info.box[0] = molecular->boxDimensions[0] * ANGSTROM_TO_NM;  // Convert Å to nm
        state.info.box[1] = molecular->boxDimensions[1] * ANGSTROM_TO_NM;  // Convert Å to nm
        state.info.box[2] = molecular->boxDimensions[2] * ANGSTROM_TO_NM;  // Convert Å to nm
        state.info.volume = state.info.box[0] * state.info.box[1] * state.info.box[2];
    } else {
        // Some minimal inputs omit CRYST1; the active box can be provided by INP instead.
        state.info.box[0] = 0.0f;
        state.info.box[1] = 0.0f;
        state.info.box[2] = 0.0f;
        state.info.volume = 0.0f;
    }

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
    // Get total number of atom types in the system
    int numTypes = state.atomTypes.atomTypes.size();
    
    state.forcefield.numTotalTypes = numTypes;
    state.forcefield.numMovementTypes = state.numMovementAtomTypes;
    
    // Resize force field arrays
    state.forcefield.ljSigma.resize(numTypes * numTypes, 0.0f);
    state.forcefield.ljEps.resize(numTypes * numTypes, 0.0f);
    
    // Convert CHARMM parameters to GROMACS format
    for (int i = 0; i < numTypes; ++i) {
        for (int j = 0; j < numTypes; ++j) {
            int idx = i * numTypes + j;
            
            // Get atom type names
            std::string type1 = state.atomTypes.getTypeName(i);
            std::string type2 = state.atomTypes.getTypeName(j);
            
            // First try to get NBFIX parameters
            auto [nbfix_params, has_nbfix] = ff.get_nbfix(type1, type2);
            
            try {
                if (has_nbfix) {
                    // Use NBFIX parameters directly
                    // Convert Rmin from Å to nm, then to sigma: sigma = Rmin / 2^(1/6)
                    const float sigma = static_cast<float>(nbfix_params.rmin / std::pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM;
                    
                    // Store parameters (eps in kJ/mol, sigma in nm)
                    state.forcefield.ljSigma[idx] = sigma;
                    state.forcefield.ljEps[idx] = static_cast<float>(nbfix_params.epsilon) * KCAL_TO_KJ;
                } else {
                    // Get LJ parameters for both types
                    auto lj1 = ff.get_lj_params(type1);
                    auto lj2 = ff.get_lj_params(type2);
                    
                    // Apply Lorentz-Berthelot combining rules
                    // sigma = (sigma1 + sigma2) / 2
                    // epsilon = sqrt(epsilon1 * epsilon2)
                    
                    // Convert CHARMM Rmin/2 to sigma: sigma = Rmin / 2^(1/6)
                    // Rmin = 2 * Rmin/2, so sigma = 2 * Rmin/2 / 2^(1/6)
                    float sigma1 = static_cast<float>(2.0 * lj1.rmin_half / std::pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM;
                    float sigma2 = static_cast<float>(2.0 * lj2.rmin_half / std::pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM;
                    state.forcefield.ljSigma[idx] = (sigma1 + sigma2) * 0.5f;
                    
                    // Convert CHARMM epsilon to kJ/mol
                    float eps1 = std::abs(lj1.epsilon) * KCAL_TO_KJ; // kcal/mol to kJ/mol
                    float eps2 = std::abs(lj2.epsilon) * KCAL_TO_KJ; // kcal/mol to kJ/mol
                    state.forcefield.ljEps[idx] = std::sqrt(eps1 * eps2);
                }
                
            } catch (const std::exception&) {
                // If parameters not found, set to zero
                state.forcefield.ljSigma[idx] = 0.0f;
                state.forcefield.ljEps[idx] = 0.0f;
            }
        }
    }
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

    if (molecular->residues.size() < numResidues || molecular->topology_residues.size() < numResidues) {
        throw std::runtime_error(
            "Inconsistent molecular data: residue count mismatch between structure and topology");
    }
    
    for (size_t i = 0; i < numResidues; ++i) {
        const auto& molRes = molecular->residues[i];
        const auto& topRes = molecular->topology_residues[i];

        if (!molRes) {
            throw std::runtime_error("Null residue encountered in molecular structure data");
        }
        
        model::MCResidue mcRes;
        mcRes.atomStart = atomStart;
        mcRes.atomCount = molRes->atom_count();
        mcRes.active = true;
        mcRes.fixed = false;
        mcRes.resname = molRes->get_resname();
        mcRes.resid = molRes->get_ires();
        
        // Initialize energy components and GCMC parameters in GROMACS units
        mcRes.energy_vdw = 0.0f;   // kJ/mole
        mcRes.energy_elec = 0.0f;  // kJ/mole
        mcRes.chemPot = 0.0f;      // kJ/mole
        mcRes.concentration = 0.0f; // mol/L
        mcRes.radius = 0.0f;       // nm
        
        mcRes.atoms.clear();
        mcRes.atoms.reserve(molRes->atom_count());

        // Convert atoms for this residue
        const auto& molAtoms = molRes->get_atoms();
        if (topRes.atoms.size() != molAtoms.size()) {
            throw std::runtime_error(
                "Inconsistent molecular data: residue atom count mismatch between structure and topology");
        }
        for (size_t j = 0; j < molAtoms.size(); j++) {
            const auto& molAtom = molAtoms[j];
            const int topAtomIdx = topRes.atoms[j];
            if (topAtomIdx < 0 || static_cast<size_t>(topAtomIdx) >= molecular->topology_atoms.size()) {
                throw std::runtime_error("Invalid topology atom index while initializing MC state");
            }
            const auto& topAtom = molecular->topology_atoms[static_cast<size_t>(topAtomIdx)];
            
            model::MCAtom mcAtom;
            // Convert coordinates from Å to nm
            mcAtom.x = molAtom->get_x() * ANGSTROM_TO_NM;
            mcAtom.y = molAtom->get_y() * ANGSTROM_TO_NM;
            mcAtom.z = molAtom->get_z() * ANGSTROM_TO_NM;
            mcAtom.charge = topAtom.charge;  // Charge unit (e) remains the same
            mcAtom.type = state.atomTypes.getOrAddType(topAtom.type);
            mcAtom.name = topAtom.name;
            mcAtom.mass = static_cast<float>(topAtom.mass);
            mcAtom.updatePosition();
            
            tempAtoms.push_back(mcAtom);
            mcRes.atoms.push_back(mcAtom);
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
