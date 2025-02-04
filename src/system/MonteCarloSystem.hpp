// src/system/MonteCarloSystem.hpp

#pragma once
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include "model/montecarlo.hpp"
#include "model/molecular.hpp"
#include "model/forcefield.hpp"

namespace pygcmc {
namespace system {

class MonteCarloSystem {
public:
    // ------------------------------------------------------------
    // Constructor / Destructor
    // ------------------------------------------------------------
    MonteCarloSystem() = default;
    ~MonteCarloSystem() = default;

    // Disable copy
    MonteCarloSystem(const MonteCarloSystem&) = delete;
    MonteCarloSystem& operator=(const MonteCarloSystem&) = delete;

    // Enable move
    MonteCarloSystem(MonteCarloSystem&&) = default;
    MonteCarloSystem& operator=(MonteCarloSystem&&) = default;

    // ------------------------------------------------------------
    // System initialization and setup
    // ------------------------------------------------------------
    void initialize(const model::MCInfo& info) {
        state.info = info;
        state.atoms.resize(info.maxAtoms);
        state.residues.resize(info.maxResidues);
    }

    void setForceField(const model::MCForceField& ff) {
        state.forcefield = ff;
    }

    void initializeForceField(const ForceField& ff);

    void addInitialResidues(const model::MCResidue* resVec, int resCount,
                           const model::MCAtom* atomVec, int atomCount);

    void initializeFromMolecular(const std::shared_ptr<model::Molecular>& molecular);

    // Movement residue management
    struct MovementMolecularInfo {
        std::shared_ptr<model::Molecular> molecular;
        int maxCopies;

        MovementMolecularInfo(std::shared_ptr<model::Molecular> mol, int max)
            : molecular(mol), maxCopies(max) {}
    };
    void addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules);

    const model::TypeMaps& getTypeMaps() const { return state.atomTypes; }

    // ------------------------------------------------------------
    // Core GCMC operations
    // ------------------------------------------------------------
    int insertResidue(const model::MCResidue& res, const model::MCAtom* atoms);
    bool removeResidue(int resIdx);
    void translateResidue(int resIdx, float dx, float dy, float dz);
    float calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const;
    float calcTotalEnergy() const;

    // ------------------------------------------------------------
    // State access
    // ------------------------------------------------------------
    const model::MCState& getState() const { return state; }
    int getActiveResidueCount() const { return state.activeResidueCount; }
    int getActiveAtomCount() const { return state.activeAtomCount; }

private:
    // ------------------------------------------------------------
    // Private helper functions
    // ------------------------------------------------------------
    void applyPBC(float& x, float& y, float& z) const;
    float getMinImageDistSqr(float dx, float dy, float dz) const;
    void updateGeometricCenter(model::MCResidue& res);

    // ------------------------------------------------------------
    // Member variables
    // ------------------------------------------------------------
    model::MCState state;
};

} // namespace system
} // namespace pygcmc

