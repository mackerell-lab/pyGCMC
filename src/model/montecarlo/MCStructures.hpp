#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>

/**
 * @file   MCStructures.hpp
 * @brief  Core data structures for GCMC simulation
 */

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief 3D Vector for positions and displacements
 */
struct Vector3 {
    double x, y, z;
    
    Vector3() : x(0), y(0), z(0) {}
    Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    Vector3 operator+(const Vector3& v) const { return Vector3(x+v.x, y+v.y, z+v.z); }
    Vector3 operator-(const Vector3& v) const { return Vector3(x-v.x, y-v.y, z-v.z); }
    Vector3 operator*(double s) const { return Vector3(x*s, y*s, z*s); }
    double dot(const Vector3& v) const { return x*v.x + y*v.y + z*v.z; }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
};

/**
 * @brief Quaternion for orientations
 */
struct Quaternion {
    double w, x, y, z;
    
    Quaternion() : w(1), x(0), y(0), z(0) {}
    Quaternion(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}
    
    void normalize() {
        double n = std::sqrt(w*w + x*x + y*y + z*z);
        if (n > 0) { w /= n; x /= n; y /= n; z /= n; }
    }
    
    // CRITICAL ADDITION: Quaternion multiplication for proper rotation composition
    Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w
        );
    }
    
    // Apply rotation to a vector
    Vector3 rotate(const Vector3& v) const {
        // Standard quaternion rotation formula: v' = q * v * q^*
        // Using the correct formula from computer graphics/robotics
        double qw = w, qx = x, qy = y, qz = z;
        double vx = v.x, vy = v.y, vz = v.z;
        
        // Rotation matrix form (verified formula)
        double qw2 = qw * qw;
        double qx2 = qx * qx;
        double qy2 = qy * qy;
        double qz2 = qz * qz;
        
        double rx = vx * (qw2 + qx2 - qy2 - qz2) + 
                   vy * 2.0 * (qx * qy - qw * qz) + 
                   vz * 2.0 * (qx * qz + qw * qy);
                   
        double ry = vx * 2.0 * (qx * qy + qw * qz) + 
                   vy * (qw2 - qx2 + qy2 - qz2) + 
                   vz * 2.0 * (qy * qz - qw * qx);
                   
        double rz = vx * 2.0 * (qx * qz - qw * qy) + 
                   vy * 2.0 * (qy * qz + qw * qx) + 
                   vz * (qw2 - qx2 - qy2 + qz2);
        
        return Vector3(rx, ry, rz);
    }
};

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
    std::vector<float> ljSigma;  // NxN matrix of sigma values (nm)
    std::vector<float> ljEps;    // NxN matrix of epsilon values (kJ/mol)
    
    // NBFIX support and mixing rules
    struct NBFixEntry {
        int type1;
        int type2;
        float sigma;  // nm
        float eps;    // kJ/mol
    };
    
    std::vector<NBFixEntry> nbfix;      // Pair-specific overrides
    std::vector<float> ljSigmaType;     // Per-type sigma values (nm)
    std::vector<float> ljEpsType;       // Per-type epsilon values (kJ/mol)
    
    enum class MixingRule { 
        None,              // Use explicit NxN matrix
        LorentzBerthelot,  // sigma_ij = (sigma_i + sigma_j)/2, eps_ij = sqrt(eps_i * eps_j)
        Geometric          // sigma_ij = sqrt(sigma_i * sigma_j), eps_ij = sqrt(eps_i * eps_j)
    };
    MixingRule mixingRule = MixingRule::None;
    
    bool ljMatrixInitialized = false;
    
    /**
     * @brief Rebuild the NxN LJ parameter matrix from per-type values and NBFIX overrides
     * 
     * This function constructs the full NxN matrix by:
     * 1. Applying mixing rules to per-type parameters
     * 2. Overriding specific pairs with NBFIX values
     */
    inline void rebuildLJMatrix() {
        const int n = numTotalTypes;
        const size_t expected = static_cast<size_t>(n) * n;
        
        // Ensure matrix has correct size
        if (ljSigma.size() != expected || ljEps.size() != expected) {
            ljSigma.assign(expected, 0.0f);
            ljEps.assign(expected, 0.0f);
        }
        
        // Build from per-type parameters with mixing rule
        if (mixingRule != MixingRule::None &&
            ljSigmaType.size() == static_cast<size_t>(n) &&
            ljEpsType.size() == static_cast<size_t>(n)) {
            
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    float sigma, eps;
                    
                    if (mixingRule == MixingRule::LorentzBerthelot) {
                        // Lorentz-Berthelot: arithmetic mean for sigma
                        sigma = 0.5f * (ljSigmaType[i] + ljSigmaType[j]);
                    } else {
                        // Geometric: geometric mean for sigma
                        sigma = std::sqrt(ljSigmaType[i] * ljSigmaType[j]);
                    }
                    
                    // Both rules use geometric mean for epsilon
                    eps = std::sqrt(ljEpsType[i] * ljEpsType[j]);
                    
                    ljSigma[i*n + j] = sigma;
                    ljEps[i*n + j] = eps;
                }
            }
        }
        
        // Apply NBFIX overrides (symmetric)
        for (const auto& p : nbfix) {
            if (p.type1 >= 0 && p.type1 < n && p.type2 >= 0 && p.type2 < n) {
                const int idx1 = p.type1 * n + p.type2;
                const int idx2 = p.type2 * n + p.type1;
                ljSigma[idx1] = ljSigma[idx2] = p.sigma;
                ljEps[idx1] = ljEps[idx2] = p.eps;
            }
        }
        
        ljMatrixInitialized = true;
    }
    
    /**
     * @brief Add an NBFIX override for a specific atom pair
     */
    inline void addNBFix(int type1, int type2, float sigma, float eps) {
        nbfix.push_back({type1, type2, sigma, eps});
        ljMatrixInitialized = false;  // Force rebuild
    }
    
    /**
     * @brief Set per-type LJ parameters
     */
    inline void setPerTypeParameters(const std::vector<float>& sigmas, 
                                     const std::vector<float>& epsilons) {
        ljSigmaType = sigmas;
        ljEpsType = epsilons;
        ljMatrixInitialized = false;  // Force rebuild
    }
};

/**
 * @brief Basic atomic properties for Monte Carlo simulation
 */
struct MCAtom {
    float x{0.0f}, y{0.0f}, z{0.0f};
    float charge{0.0f};
    float mass{1.0f};  // Atom mass in amu
    int   type{-1};
    std::string name;  // Atom name (e.g., "O", "H1", "H2")
    Vector3 position;  // Position as Vector3 (for convenience)
    
    // Helper to update position from x,y,z
    void updatePosition() {
        position = Vector3(x, y, z);
    }
    
    // Helper to set x,y,z from position
    void setFromPosition() {
        x = static_cast<float>(position.x);
        y = static_cast<float>(position.y);
        z = static_cast<float>(position.z);
    }
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
    std::string resname;  // Residue name (e.g., "TIP3", "WAT", etc.)
    int resid{-1};        // Residue ID
    std::vector<MCAtom> atoms;  // Atoms in this residue
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