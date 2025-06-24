#pragma once

#ifndef PYGCMC_MODEL_MONTECARLO_STRUCTURES_HPP
#define PYGCMC_MODEL_MONTECARLO_STRUCTURES_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <array>
#include <stdexcept>
#include <cmath>

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

    size_t size() const { return atomTypes.size(); }

    void clear() {
        atomTypes.clear();
        atomTypeIndices.clear();
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

    // Reserved max capacity
    int    maxResidues;
    int    maxAtoms;
    int    maxTypes;

    // Global parameters
    float  volume;
    uint64_t seed;

    // CHARMM-style switching function parameters
    bool use_switching{false};
    float r_on{1.0f};
    float r_off{1.2f};

    /**
     * @brief Statistics for Monte Carlo moves
     */
    struct Statistics {
        int totalMoves{0};
        int acceptedMoves{0};
        int insertionAttempts{0};
        int acceptedInsertions{0};
        int deletionAttempts{0};
        int acceptedDeletions{0};

        double getAcceptanceRate() const {
            return (totalMoves > 0) ? static_cast<double>(acceptedMoves) / totalMoves : 0.0;
        }

        double getInsertionRate() const {
            return (insertionAttempts > 0) ? static_cast<double>(acceptedInsertions) / insertionAttempts : 0.0;
        }

        double getDeletionRate() const {
            return (deletionAttempts > 0) ? static_cast<double>(acceptedDeletions) / deletionAttempts : 0.0;
        }

        void reset() {
            totalMoves = 0;
            acceptedMoves = 0;
            insertionAttempts = 0;
            acceptedInsertions = 0;
            deletionAttempts = 0;
            acceptedDeletions = 0;
        }
    } stats;

    // Physical constants
    static constexpr float BOLTZMANN = 0.00831446f;
    static constexpr float MOLES_TO_MOLECULES = 0.0006023f;
    static constexpr float MOLECULES_TO_MOLES = 1660.539f;

    // Temperature management methods
    void setTemperature(float temperature) {
        if (temperature <= 0.0f) {
            throw std::invalid_argument("Temperature must be positive");
        }
        beta = 1.0f / (BOLTZMANN * temperature);
    }

    float getTemperature() const {
        return (beta > 0.0f) ? 1.0f / (BOLTZMANN * beta) : 0.0f;
    }
};

/**
 * @brief Force field parameters for Monte Carlo simulation
 */
struct MCForceField {
    int numTotalTypes;
    int numMovementTypes;

    std::vector<float> ljSigma;
    std::vector<float> ljEps;
};

/**
 * @brief Basic atomic properties for Monte Carlo simulation
 */
struct MCAtom {
    float x, y, z;
    float charge;
    int   type;

    MCAtom() : x(0.0f), y(0.0f), z(0.0f), charge(0.0f), type(-1) {}
    MCAtom(float x_, float y_, float z_, float charge_, int type_) 
        : x(x_), y(y_), z(z_), charge(charge_), type(type_) {}

    void setPosition(float x_, float y_, float z_) {
        x = x_;
        y = y_;
        z = z_;
    }

    std::array<float, 3> getPosition() const {
        return {x, y, z};
    }
};

/**
 * @brief Molecular unit for GCMC simulation
 */
struct MCResidue {
    // Basic properties
    int   atomStart;
    int   atomCount;
    bool  active;
    bool  fixed;
    float center[3];

    // Energy components
    float energy_vdw;
    float energy_elec;
        
    // GCMC parameters
    float concentration;
    float chemPot;
    int   type;
    float radius;

    MCResidue() : atomStart(-1), atomCount(0), active(false), fixed(false),
                  center{0.0f, 0.0f, 0.0f}, energy_vdw(0.0f), energy_elec(0.0f),
                  concentration(0.0f), chemPot(0.0f), type(-1), radius(0.0f) {}

    float getTotalEnergy() const {
        return energy_vdw + energy_elec;
    }

    void setCenter(float x, float y, float z) {
        center[0] = x;
        center[1] = y;
        center[2] = z;
    }

    std::array<float, 3> getCenter() const {
        return {center[0], center[1], center[2]};
    }

    bool isValid() const {
        return atomStart >= 0 && atomCount > 0 && type >= 0;
    }
};

/**
 * @brief Information about movement residues in the system
 */
struct MCMovementResidueInfo {
    int startIndex;
    int activeCount;
    int totalCount;
    std::string resName;

    MCMovementResidueInfo() : startIndex(-1), activeCount(0), totalCount(0) {}
    MCMovementResidueInfo(int start, int active, int total, const std::string& name)
        : startIndex(start), activeCount(active), totalCount(total), resName(name) {}

    double getActiveFraction() const {
        return (totalCount > 0) ? static_cast<double>(activeCount) / totalCount : 0.0;
    }
};

/**
 * @brief Ewald energy components
 */
struct EwaldEnergy {
    double real_space{0.0};
    double reciprocal{0.0};
    double self{0.0};
    double total{0.0};

    void reset() {
        real_space = 0.0;
        reciprocal = 0.0;
        self = 0.0;
        total = 0.0;
    }

    void updateTotal() {
        total = real_space + reciprocal + self;
    }
};

} // namespace montecarlo
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MONTECARLO_STRUCTURES_HPP