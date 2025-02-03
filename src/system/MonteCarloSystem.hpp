// src/system/MonteCarloSystem.hpp

#pragma once
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include "model/montecarlo.hpp"
#include "model/molecular.hpp"

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

    const model::TypeMaps& getTypeMaps() const { return typeMaps; }

    // ------------------------------------------------------------
    // Core GCMC operations
    // ------------------------------------------------------------
    int insertResidue(const model::MCResidue& res, const model::MCAtom* atoms);
    bool removeResidue(int resIdx);
    void translateResidue(int resIdx, float dx, float dy, float dz);

    // ------------------------------------------------------------
    // Energy calculation
    // ------------------------------------------------------------
    float calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const;
    float calcTotalEnergy() const;

    // ------------------------------------------------------------
    // Periodic boundary conditions
    // ------------------------------------------------------------
    void applyPBC(float& x, float& y, float& z) const;
    float getMinImageDistSqr(float dx, float dy, float dz) const;

    // ------------------------------------------------------------
    // System state access
    // ------------------------------------------------------------
    const model::MCState& getState() const { return state; }
    model::MCState& getState() { return state; }

    int getActiveAtomCount() const { return state.activeAtomCount; }
    int getActiveResidueCount() const { return state.activeResidueCount; }

private:
    // Helper methods
    void updateGeometricCenter(model::MCResidue& res);

    // System state
    model::MCState state;
    model::TypeMaps typeMaps;
};

} // namespace system
} // namespace pygcmc

