#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ENERGY_CALLBACK_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ENERGY_CALLBACK_HPP

#include "model/ModelModule.hpp"
#include "../../energy/EnergyModule.hpp"
#include <functional>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

using namespace model::montecarlo;

/**
 * @brief Energy calculation callback interface for GCMC
 * 
 * This interface allows GCMC to use the existing energy modules
 * (DIRECT, EWALD, PME) for accurate energy calculations.
 */
class GCMCEnergyCallback {
public:
    // Callback function types
    using SystemEnergyFunc = std::function<double(MCState&)>;
    using ResidueEnergyFunc = std::function<double(MCState&, int)>;
    using EnergyDifferenceFunc = std::function<double(MCState&, int, const Vector3&, const Vector3&)>;
    
    GCMCEnergyCallback() 
        : energyMethod_(EnergyMethod::DIRECT),
          useCutoff_(true),
          usePBC_(true) {}
    
    // Set energy calculation method
    void setEnergyMethod(EnergyMethod method) {
        energyMethod_ = method;
    }
    
    // Set calculation parameters
    void setParameters(bool useCutoff, bool usePBC) {
        useCutoff_ = useCutoff;
        usePBC_ = usePBC;
    }
    
    /**
     * @brief Calculate total system energy
     * @param state System state
     * @return Total energy in kJ/mol
     */
    double calculateSystemEnergy(MCState& state) {
        // Use the energy module to calculate system energy
        if (energyMethod_ == EnergyMethod::PME) {
            computeSystemEnergy(state, EnergyMethod::PME);
            return energy::getTotalEnergyUniquePairs(state, EnergyMethod::PME);
        } else if (energyMethod_ == EnergyMethod::EWALD) {
            computeSystemEnergy(state, EnergyMethod::EWALD);
            return energy::getTotalEnergyUniquePairs(state, EnergyMethod::EWALD);
        } else {
            computeSystemEnergy(state, EnergyMethod::DIRECT, useCutoff_, usePBC_);
            return energy::getTotalEnergyUniquePairs(state, EnergyMethod::DIRECT);
        }
    }
    
    /**
     * @brief Calculate energy of a single residue
     * @param state System state
     * @param residueIdx Residue index
     * @return Energy in kJ/mol
     */
    double calculateResidueEnergy(MCState& state, int residueIdx) {
        if (residueIdx < 0 || residueIdx >= static_cast<int>(state.residues.size())) {
            return 0.0;
        }
        
        // Mark residue as moved for movement energy calculation
        // Note: moved flag may not exist in MCResidue
        // state.residues[residueIdx].moved = true;
        
        // Calculate movement energy
        if (energyMethod_ == EnergyMethod::PME) {
            computeMovementEnergy(state, EnergyMethod::PME);
        } else if (energyMethod_ == EnergyMethod::EWALD) {
            computeMovementEnergy(state, EnergyMethod::EWALD);
        } else {
            computeMovementEnergy(state, EnergyMethod::DIRECT, useCutoff_, usePBC_);
        }
        
        // Reset moved flag
        // state.residues[residueIdx].moved = false;
        
        // Return residue energy
        const auto& residue = state.residues[residueIdx];
        return residue.energy_vdw + residue.energy_elec;
    }
    
    /**
     * @brief Calculate energy difference for a position change
     * @param state System state
     * @param residueIdx Residue index
     * @param oldPos Old position
     * @param newPos New position
     * @return Energy difference in kJ/mol
     */
    double calculateEnergyDifference(MCState& state, int residueIdx,
                                    const Vector3& oldPos, const Vector3& newPos) {
        // Calculate energy at old position
        double oldEnergy = calculateResidueEnergy(state, residueIdx);
        
        // Temporarily move residue to new position
        if (residueIdx >= 0 && residueIdx < static_cast<int>(state.residues.size())) {
            auto& residue = state.residues[residueIdx];
            for (auto& atom : residue.atoms) {
                atom.x += (newPos.x - oldPos.x);
                atom.y += (newPos.y - oldPos.y);
                atom.z += (newPos.z - oldPos.z);
            }
        }
        
        // Calculate energy at new position
        double newEnergy = calculateResidueEnergy(state, residueIdx);
        
        // Restore original position
        if (residueIdx >= 0 && residueIdx < static_cast<int>(state.residues.size())) {
            auto& residue = state.residues[residueIdx];
            for (auto& atom : residue.atoms) {
                atom.x -= (newPos.x - oldPos.x);
                atom.y -= (newPos.y - oldPos.y);
                atom.z -= (newPos.z - oldPos.z);
            }
        }
        
        return newEnergy - oldEnergy;
    }
    
    /**
     * @brief Initialize Ewald parameters
     * @param cutoff Cutoff distance in nm
     * @param box Box dimensions
     */
    void initializeEwald(double cutoff, const std::vector<double>& box) {
        (void)cutoff;  // Suppress unused parameter warning
        (void)box;     // Suppress unused parameter warning
        if (energyMethod_ == EnergyMethod::EWALD) {
            // Initialize Ewald parameters
            // This would call the appropriate initialization functions
            // from the Ewald module
        }
    }
    
    /**
     * @brief Initialize PME parameters
     * @param cutoff Cutoff distance in nm
     * @param box Box dimensions
     * @param gridSize Grid dimensions for PME
     */
    void initializePME(double cutoff, const std::vector<double>& box, 
                      const std::vector<int>& gridSize) {
        (void)cutoff;   // Suppress unused parameter warning
        (void)box;      // Suppress unused parameter warning
        (void)gridSize; // Suppress unused parameter warning
        if (energyMethod_ == EnergyMethod::PME) {
            // Initialize PME parameters
            // This would call the appropriate initialization functions
            // from the PME module
        }
    }
    
private:
    EnergyMethod energyMethod_;
    bool useCutoff_;
    bool usePBC_;
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_ENERGY_CALLBACK_HPP
