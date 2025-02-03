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
    // Type mapping system
    struct TypeMaps {
        std::vector<std::string> atomTypes;  // Index -> Type string mapping
        std::unordered_map<std::string, int> atomTypeIndices;  // Type string -> Index mapping
        
        int getOrAddType(const std::string& type) {
            auto it = atomTypeIndices.find(type);
            if (it != atomTypeIndices.end()) {
                return it->second;
            }
            int newIndex = atomTypes.size();
            atomTypes.push_back(type);
            atomTypeIndices[type] = newIndex;
            return newIndex;
        }
        
        std::string getTypeName(int index) const {
            if (index >= 0 && static_cast<size_t>(index) < atomTypes.size()) {
                return atomTypes[index];
            }
            return "";
        }
    };

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
    void initialize(const GCMCInfo& info) {
        state.info = info;
        state.atoms.resize(info.maxAtoms);
        state.residues.resize(info.maxResidues);
    }

    void setForceField(const ForceField& ff) {
        state.forcefield = ff;
    }

    void addInitialResidues(const Residue* resVec, int resCount,
                           const Atom* atomVec, int atomCount);

    void initializeFromMolecular(const std::shared_ptr<pygcmc::model::Molecular>& molecular);

    const TypeMaps& getTypeMaps() const { return typeMaps; }

    // ------------------------------------------------------------
    // Core GCMC operations
    // ------------------------------------------------------------
    int insertResidue(const Residue& res, const Atom* atoms);
    bool removeResidue(int resIdx);
    void translateResidue(int resIdx, float dx, float dy, float dz);

    // ------------------------------------------------------------
    // Energy calculation
    // ------------------------------------------------------------
    float calcNonBondedEnergy(const Residue& res1, const Residue& res2) const;
    float calcTotalEnergy() const;

    // ------------------------------------------------------------
    // Periodic boundary conditions
    // ------------------------------------------------------------
    void applyPBC(float& x, float& y, float& z) const;
    float getMinImageDistSqr(float dx, float dy, float dz) const;

    // ------------------------------------------------------------
    // System state access
    // ------------------------------------------------------------
    const SystemState& getState() const { return state; }
    SystemState& getState() { return state; }

    int getActiveAtomCount() const { return state.activeAtomCount; }
    int getActiveResidueCount() const { return state.activeResidueCount; }

private:
    // Helper methods
    void updateCenterOfMass(Residue& res);

    // System state
    SystemState state;
    TypeMaps typeMaps;
};

} // namespace gcmc

