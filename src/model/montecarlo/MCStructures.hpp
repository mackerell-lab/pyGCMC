#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

/**
 * @file   MCStructures.hpp
 * @brief  Core data structures for GCMC simulation
 *
 * Pure data structures with minimal methods. Units:
 * - Length: nanometers (nm), Energy: kilojoules per mole (kJ/mol)
 * - Charge: electron charge (e), Temperature: Kelvin (K)
 */

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Type mapping system for atom and residue types
 */
struct TypeMaps {
    std::vector<std::string> atomTypes;
    std::unordered_map<std::string, int> atomTypeIndices;
    
    int getOrAddType(const std::string& type) {
        auto it = atomTypeIndices.find(type);
        if (it != atomTypeIndices.end()) {
            return it->second;
        }
        int newIndex = static_cast<int>(atomTypes.size());
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

/**
 * @brief Global Monte Carlo simulation parameters
 */
struct MCInfo {
    int    mcSteps{0};
    float  box[3]{-1.0f};
    float  cutoff{1.5f};
    float  beta{0.0f};
    int    maxResidues{0};
    int    maxAtoms{0};
    int    maxTypes{0};
    float  volume{0.0f};
    uint64_t seed{0};
    bool use_switching{false};
    float r_on{1.0f};
    float r_off{1.2f};

    struct Statistics {
        int totalMoves{0};
        int acceptedMoves{0};
        int insertionAttempts{0};
        int acceptedInsertions{0};
        int deletionAttempts{0};
        int acceptedDeletions{0};
    } stats;

    static constexpr float BOLTZMANN = 0.00831446f;
    static constexpr float MOLES_TO_MOLECULES = 0.0006023f;
    static constexpr float MOLECULES_TO_MOLES = 1660.539f;

    void setTemperature(float temperature) {
        beta = 1.0f / (BOLTZMANN * temperature);
    }
};

/**
 * @brief Force field parameters for Monte Carlo simulation
 */
struct MCForceField {
    int numTotalTypes{0};
    int numMovementTypes{0};
    std::vector<float> ljSigma;
    std::vector<float> ljEps;
};

/**
 * @brief Basic atomic properties for Monte Carlo simulation
 */
struct MCAtom {
    float x{0.0f}, y{0.0f}, z{0.0f};
    float charge{0.0f};
    int   type{-1};
};

/**
 * @brief Molecular unit for GCMC simulation
 */
struct MCResidue {
    int   atomStart{-1};
    int   atomCount{0};
    bool  active{false};
    bool  fixed{false};
    float center[3]{0.0f};
    float energy_vdw{0.0f};
    float energy_elec{0.0f};
    float concentration{0.0f};
    float chemPot{0.0f};
    int   type{-1};
    float radius{0.0f};
};

/**
 * @brief Information about movement residues in the system
 */
struct MCMovementResidueInfo {
    int startIndex{-1};
    int activeCount{0};
    int totalCount{0};
    std::string resName;
};

/**
 * @brief Ewald energy components
 */
struct EwaldEnergy {
    double real_space{0.0};
    double reciprocal{0.0};
    double self{0.0};
    double total{0.0};
    
    void reset() { real_space = reciprocal = self = total = 0.0; }
    void updateTotal() { total = real_space + reciprocal + self; }
};

} // namespace montecarlo
} // namespace model
} // namespace pygcmc