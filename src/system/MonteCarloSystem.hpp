// src/system/MonteCarloSystem.hpp

#pragma once
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include "model/montecarlo.hpp"
#include "model/molecular.hpp"

namespace gcmc {

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
    void initialize(const MCInfo& info) {
        state.info = info;
        state.atoms.resize(info.maxAtoms);
        state.residues.resize(info.maxResidues);
    }

    void setForceField(const MCForceField& ff) {
        state.forcefield = ff;
    }

    void addInitialResidues(const MCResidue* resVec, int resCount,
                           const MCAtom* atomVec, int atomCount);

    void initializeFromMolecular(const std::shared_ptr<pygcmc::model::Molecular>& molecular);

    // Movement residue management
    struct MovementMolecularInfo {
        std::shared_ptr<pygcmc::model::Molecular> molecular;
        int maxCopies;
    };
    void addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules);

    const TypeMaps& getTypeMaps() const { return typeMaps; }

    // ------------------------------------------------------------
    // Core GCMC operations
    // ------------------------------------------------------------
    int insertResidue(const MCResidue& res, const MCAtom* atoms);
    bool removeResidue(int resIdx);
    void translateResidue(int resIdx, float dx, float dy, float dz);

    // ------------------------------------------------------------
    // Energy calculation
    // ------------------------------------------------------------
    float calcNonBondedEnergy(const MCResidue& res1, const MCResidue& res2) const;
    float calcTotalEnergy() const;

    // ------------------------------------------------------------
    // Periodic boundary conditions
    // ------------------------------------------------------------
    void applyPBC(float& x, float& y, float& z) const;
    float getMinImageDistSqr(float dx, float dy, float dz) const;

    // ------------------------------------------------------------
    // System state access
    // ------------------------------------------------------------
    const MCState& getState() const { return state; }
    MCState& getState() { return state; }

    int getActiveAtomCount() const { return state.activeAtomCount; }
    int getActiveResidueCount() const { return state.activeResidueCount; }

private:
    // Helper methods
    void updateGeometricCenter(MCResidue& res);

    // System state
    MCState state;
    TypeMaps typeMaps;
};

} // namespace gcmc

