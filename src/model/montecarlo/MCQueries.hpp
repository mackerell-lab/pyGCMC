#pragma once

#ifndef PYGCMC_MODEL_MONTECARLO_QUERIES_HPP
#define PYGCMC_MODEL_MONTECARLO_QUERIES_HPP

#include "MCStructures.hpp"
#include <string>
#include <sstream>
#include <optional>

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Static query and validation methods for Monte Carlo state
 */
class MCQueries {
public:
    // === MCInfo Queries ===
    static bool isValidTemperature(const MCInfo& info) {
        return info.beta > 0.0f;
    }

    static bool isValidBox(const MCInfo& info) {
        return info.box[0] > 0.0f && info.box[1] > 0.0f && info.box[2] > 0.0f;
    }

    static float getTemperature(const MCInfo& info) {
        return (info.beta > 0.0f) ? 1.0f / (MCInfo::BOLTZMANN * info.beta) : 0.0f;
    }

    static float getVolume(const MCInfo& info) {
        return info.box[0] * info.box[1] * info.box[2];
    }

    static double getAcceptanceRate(const MCInfo::Statistics& stats) {
        return (stats.totalMoves > 0) ? static_cast<double>(stats.acceptedMoves) / stats.totalMoves : 0.0;
    }

    // === MCForceField Queries ===
    static bool isValidForceField(const MCForceField& ff) {
        if (ff.numTotalTypes <= 0) return false;
        size_t expected_size = ff.numTotalTypes * ff.numTotalTypes;
        return ff.ljSigma.size() == expected_size && ff.ljEps.size() == expected_size;
    }

    static std::pair<float, float> getLJParams(const MCForceField& ff, int type1, int type2) {
        if (type1 >= ff.numTotalTypes || type2 >= ff.numTotalTypes || type1 < 0 || type2 < 0) {
            return std::make_pair(0.0f, 0.0f);
        }
        int index = type1 * ff.numTotalTypes + type2;
        return std::make_pair(ff.ljSigma[index], ff.ljEps[index]);
    }

    static bool hasLJParams(const MCForceField& ff, int type1, int type2) {
        if (type1 >= ff.numTotalTypes || type2 >= ff.numTotalTypes || type1 < 0 || type2 < 0) {
            return false;
        }
        int index = type1 * ff.numTotalTypes + type2;
        return ff.ljSigma[index] > 0.0f || ff.ljEps[index] > 0.0f;
    }

    // === Atom Queries ===
    static bool isValidAtom(const MCAtom& atom) {
        return atom.type >= 0;
    }

    static float getDistance(const MCAtom& atom1, const MCAtom& atom2) {
        float dx = atom1.x - atom2.x;
        float dy = atom1.y - atom2.y;
        float dz = atom1.z - atom2.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    static std::optional<int> findAtomByType(const std::vector<MCAtom>& atoms, int activeCount, int type) {
        for (int i = 0; i < activeCount && i < static_cast<int>(atoms.size()); ++i) {
            if (atoms[i].type == type) {
                return i;
            }
        }
        return std::nullopt;
    }

    static int countAtomsByType(const std::vector<MCAtom>& atoms, int activeCount, int type) {
        int count = 0;
        for (int i = 0; i < activeCount && i < static_cast<int>(atoms.size()); ++i) {
            if (atoms[i].type == type) {
                count++;
            }
        }
        return count;
    }

    // === Residue Queries ===
    static bool isValidResidue(const MCResidue& residue) {
        return residue.atomStart >= 0 && residue.atomCount > 0 && residue.type >= 0;
    }

    static bool isActiveResidue(const MCResidue& residue) {
        return residue.active && isValidResidue(residue);
    }

    static std::optional<int> findResidueByType(const std::vector<MCResidue>& residues, int activeCount, int type) {
        for (int i = 0; i < activeCount && i < static_cast<int>(residues.size()); ++i) {
            if (residues[i].type == type && residues[i].active) {
                return i;
            }
        }
        return std::nullopt;
    }

    static int countActiveResidues(const std::vector<MCResidue>& residues, int activeCount) {
        int count = 0;
        for (int i = 0; i < activeCount && i < static_cast<int>(residues.size()); ++i) {
            if (residues[i].active) {
                count++;
            }
        }
        return count;
    }

    static int countResiduesByType(const std::vector<MCResidue>& residues, int activeCount, int type) {
        int count = 0;
        for (int i = 0; i < activeCount && i < static_cast<int>(residues.size()); ++i) {
            if (residues[i].type == type && residues[i].active) {
                count++;
            }
        }
        return count;
    }

    static float getTotalResidueEnergy(const std::vector<MCResidue>& residues, int activeCount) {
        float total = 0.0f;
        for (int i = 0; i < activeCount && i < static_cast<int>(residues.size()); ++i) {
            if (residues[i].active) {
                total += residues[i].getTotalEnergy();
            }
        }
        return total;
    }

    // === TypeMaps Queries ===
    static bool hasType(const TypeMaps& typeMaps, const std::string& type) {
        return typeMaps.atomTypeIndices.find(type) != typeMaps.atomTypeIndices.end();
    }

    static std::optional<int> getTypeIndex(const TypeMaps& typeMaps, const std::string& type) {
        auto it = typeMaps.atomTypeIndices.find(type);
        if (it != typeMaps.atomTypeIndices.end()) {
            return it->second;
        }
        return std::nullopt;
    }

    static bool isValidTypeIndex(const TypeMaps& typeMaps, int index) {
        return index >= 0 && static_cast<size_t>(index) < typeMaps.atomTypes.size();
    }

    // === Movement Residue Queries ===
    static std::optional<size_t> findMovementResidueByName(const std::vector<MCMovementResidueInfo>& movementResidues,
                                                         const std::string& resName) {
        for (size_t i = 0; i < movementResidues.size(); ++i) {
            if (movementResidues[i].resName == resName) {
                return i;
            }
        }
        return std::nullopt;
    }

    static int getTotalMovementResidues(const std::vector<MCMovementResidueInfo>& movementResidues) {
        int total = 0;
        for (const auto& info : movementResidues) {
            total += info.totalCount;
        }
        return total;
    }

    static int getActiveMovementResidues(const std::vector<MCMovementResidueInfo>& movementResidues) {
        int active = 0;
        for (const auto& info : movementResidues) {
            active += info.activeCount;
        }
        return active;
    }

    // === System Validation ===
    static bool isConsistentState(const std::vector<MCAtom>& atoms, const std::vector<MCResidue>& residues,
                                int activeAtomCount, int activeResidueCount) {
        // Check basic bounds
        if (activeAtomCount < 0 || activeResidueCount < 0) return false;
        if (activeAtomCount > static_cast<int>(atoms.size())) return false;
        if (activeResidueCount > static_cast<int>(residues.size())) return false;
        
        // Check residue-atom consistency
        for (int i = 0; i < activeResidueCount; ++i) {
            const auto& residue = residues[i];
            if (!residue.active) continue;
            if (residue.atomStart < 0 || residue.atomCount <= 0) return false;
            if (residue.atomStart + residue.atomCount > activeAtomCount) return false;
        }
        
        return true;
    }

    // === String Representations ===
    static std::string toString(const MCInfo& info) {
        std::stringstream ss;
        ss << "MCInfo: T=" << getTemperature(info) << "K, ";
        ss << "Box=[" << info.box[0] << "," << info.box[1] << "," << info.box[2] << "], ";
        ss << "Cutoff=" << info.cutoff << "nm";
        return ss.str();
    }

    static std::string toString(const MCAtom& atom) {
        std::stringstream ss;
        ss << "Atom: pos=(" << atom.x << "," << atom.y << "," << atom.z << "), ";
        ss << "q=" << atom.charge << ", type=" << atom.type;
        return ss.str();
    }

    static std::string toString(const MCResidue& residue) {
        std::stringstream ss;
        ss << "Residue: atoms=[" << residue.atomStart << ":" << residue.atomCount << "], ";
        ss << "active=" << (residue.active ? "true" : "false") << ", ";
        ss << "energy=" << residue.getTotalEnergy() << "kJ/mol, type=" << residue.type;
        return ss.str();
    }

    static std::string toString(const EwaldEnergy& ewald) {
        std::stringstream ss;
        ss << "Ewald: real=" << ewald.real_space << ", recip=" << ewald.reciprocal;
        ss << ", self=" << ewald.self << ", total=" << ewald.total << " kJ/mol";
        return ss.str();
    }
};

} // namespace montecarlo
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MONTECARLO_QUERIES_HPP