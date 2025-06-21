#include "MCComposite.hpp"
#include "MCMovementTypeCollector.hpp"
#include "MCMovementReindexer.hpp"
#include "MCMovementBuilder.hpp"

namespace pygcmc {
namespace system {
namespace montecarlo {

MCComposite::MCComposite() {
    // Initialize all subsystems
}

void MCComposite::initialize(const model::MCInfo& info) {
    state_.info = info;
    state_.atoms.resize(info.maxAtoms);
    state_.residues.resize(info.maxResidues);
}

void MCComposite::setForceField(const model::MCForceField& ff) {
    state_.forcefield = ff;
}

void MCComposite::initializeForceField(const model::ForceField& ff) {
    initializer_.initializeForceField(state_, ff);
}

void MCComposite::initializeFromMolecular(const std::shared_ptr<model::Molecular>& molecular) {
    // Store molecular system for later parameter validation
    molecular_ = molecular;
    
    // Initialize state from molecular system
    initializer_.initializeFromMolecular(state_, molecular);
}

void MCComposite::addInitialResidues(const model::MCResidue* resVec, int resCount,
                                    const model::MCAtom* atomVec, int atomCount) {
    core_.addInitialResidues(state_, resVec, resCount, atomVec, atomCount);
}

void MCComposite::addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules) {
    // Step 1: Collect types from new molecules
    MCMovementTypeCollector typeCollector;
    model::TypeMaps preResidueTypes = state_.residueTypes;
    model::TypeMaps preAtomTypes = state_.atomTypes;
    
    std::vector<int> newMovementAtomTypes = typeCollector.collectMovementTypes(
        molecules, preResidueTypes, preAtomTypes);

    // Step 2: Create a new MCState using the new type maps
    model::MCState newState;
    newState.info = state_.info;
    newState.forcefield = state_.forcefield;
    newState.residueTypes = preResidueTypes;
    newState.atomTypes = preAtomTypes;
    
    // Preserve previous movement residues information
    newState.movementResidues = state_.movementResidues;
    newState.movementAtomTypes = state_.movementAtomTypes;
    newState.numMovementAtomTypes = state_.numMovementAtomTypes;

    // Add new movement atom types
    for (int typeIdx : newMovementAtomTypes) {
        if (std::find(newState.movementAtomTypes.begin(), newState.movementAtomTypes.end(), typeIdx) 
            == newState.movementAtomTypes.end()) {
            newState.movementAtomTypes.push_back(typeIdx);
        }
    }
    newState.numMovementAtomTypes = newState.movementAtomTypes.size();

    // Step 3: Reindex existing residues and atoms
    MCMovementReindexer reindexer;
    auto [reindexedResidues, reindexedAtoms] = reindexer.reindexExistingResidues(
        state_, newState.residueTypes, newState.atomTypes);

    // Step 4: Build final state
    MCMovementBuilder builder;
    builder.buildFinalState(newState, molecules, typeCollector.getProcessedResidueNames(),
                           reindexedResidues, reindexedAtoms);

    // Update the class member state with the new state
    state_ = std::move(newState);
}

int MCComposite::insertResidue(const model::MCResidue& res, const model::MCAtom* atoms) {
    return core_.insertResidue(state_, res, atoms);
}

bool MCComposite::removeResidue(int resIdx) {
    return core_.removeResidue(state_, resIdx);
}

void MCComposite::translateResidue(int resIdx, float dx, float dy, float dz) {
    core_.translateResidue(state_, resIdx, dx, dy, dz);
}

float MCComposite::calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const {
    (void)res1;  // Suppress unused parameter warning
    (void)res2;  // Suppress unused parameter warning
    float energy = 0.0f;
    // TODO: Implement LJ + Coulomb with periodic boundary conditions
    return energy;
}

float MCComposite::calcTotalEnergy() const {
    float energy = 0.0f;
    // TODO: Implement total energy calculation
    return energy;
}

void MCComposite::setSwitchingFunction(bool enable, float r_on, float r_off) {
    switching_.setSwitchingFunction(state_, enable, r_on, r_off);
}

float MCComposite::calculateSwitchingFunction(float r) const {
    return switching_.calculateSwitchingFunction(state_, r);
}

void MCComposite::validateParameters(const model::ForceField& ff, const std::shared_ptr<model::Molecular>& molecular) {
    initializer_.validateParameters(ff, molecular);
}

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 