// src/data/residue.hpp

#ifndef PYGCMC_DATA_RESIDUE_HPP
#define PYGCMC_DATA_RESIDUE_HPP

#include <vector>
#include <memory>
#include <string>
#include <array>
#include <algorithm>
#include <unordered_map>
#include "atom.hpp"

namespace pygcmc {
namespace data {

/**
 * @brief Residue class following CHARMM naming conventions
 */
class Residue {
public:
    // Secondary structure types
    enum class SecondaryStructure {
        NONE,
        HELIX_ALPHA_RIGHT,    // Right-handed alpha (default)
        HELIX_OMEGA_RIGHT,    // Right-handed omega
        HELIX_PI_RIGHT,       // Right-handed pi
        HELIX_GAMMA_RIGHT,    // Right-handed gamma
        HELIX_310_RIGHT,      // Right-handed 3/10
        HELIX_ALPHA_LEFT,     // Left-handed alpha
        HELIX_OMEGA_LEFT,     // Left-handed omega
        HELIX_GAMMA_LEFT,     // Left-handed gamma
        HELIX_27_RIBBON,      // 2/7 ribbon/helix
        HELIX_POLYPROLINE,    // Polyproline
        SHEET_PARALLEL,       // Parallel beta sheet
        SHEET_ANTIPARALLEL    // Anti-parallel beta sheet
    };

    // Sheet strand information
    struct SheetStrand {
        int strandNumber;     // Strand number in current sheet
        std::string sheetId;  // Sheet identifier
        int numStrands;       // Number of strands in current sheet
        int sense;            // 0: first, 1: parallel, -1: antiparallel
        // Hydrogen bond information
        std::string atom1;    // First atom name
        std::string atom2;    // Second atom name
        int resnum1;          // First residue number
        int resnum2;          // Second residue number
    };

    // Disulfide bond information
    struct SSBond {
        int serialNumber;     // Bond serial number
        std::string partner;  // Partner residue identifier
        char partnerChain;    // Partner chain identifier
        int partnerResnum;    // Partner residue number
        char partnerInscode;  // Partner insertion code
        double bondLength;    // Length of the disulfide bond
    };

    // Constructors
    Residue() : 
        resname(""),      // Residue name (e.g., ALA)
        ires(0),         // Residue number (IRES)
        segid(""),       // Segment ID (SEGID)
        iseg(0),         // Segment number (ISEG)
        igro(0),         // Group number (IGRO)
        chain(' '),      // Chain ID
        inscode(' '),    // Insertion code
        move(1),         // Movement flag (MOVE)
        ignore(0),       // Ignore flag (IGNORE)
        constrain(0)     // Constraint flag (CONSTRAIN)
    {}

    Residue(const std::string& resname, int ires, 
            const std::string& segid = "", int iseg = 0,
            char chain = ' ', char inscode = ' ') :
        resname(resname),
        ires(ires),
        segid(segid),
        iseg(iseg),
        igro(0),
        chain(chain),
        inscode(inscode),
        move(1),
        ignore(0),
        constrain(0)
    {}

    // CHARMM standard getters
    const std::string& getResname() const noexcept { return resname; }
    int getIres() const noexcept { return ires; }
    const std::string& getSegid() const noexcept { return segid; }
    int getIseg() const noexcept { return iseg; }
    char getChain() const noexcept { return chain; }
    char getInscode() const noexcept { return inscode; }
    
    // Atom management
    void addAtom(const Atom& atom) {
        // Verify atom belongs to this residue
        if (atom.getResname() != resname || atom.getIres() != ires ||
            atom.getSegid() != segid || atom.getIseg() != iseg) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        atoms.push_back(std::make_shared<Atom>(atom));
    }

    void addAtom(std::shared_ptr<Atom> atom) {
        if (!atom) return;
        if (atom->getResname() != resname || atom->getIres() != ires ||
            atom->getSegid() != segid || atom->getIseg() != iseg) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        atoms.push_back(atom);
    }

    const std::vector<std::shared_ptr<Atom>>& getAtoms() const noexcept { 
        return atoms; 
    }

    std::shared_ptr<Atom> findAtom(const std::string& type) const {
        auto it = std::find_if(atoms.begin(), atoms.end(),
            [&type](const std::shared_ptr<Atom>& atom) {
                return atom && atom->getType() == type;
            });
        return (it != atoms.end()) ? *it : nullptr;
    }

    // Utility methods
    size_t atomCount() const noexcept {
        return atoms.size();
    }

    bool isValid() const {
        return !resname.empty() && ires > 0 && !segid.empty() &&
               std::all_of(atoms.begin(), atoms.end(),
                          [](const std::shared_ptr<Atom>& atom) {
                              return atom && atom->isValid();
                          });
    }

    std::array<double, 3> centerOfMass() const {
        std::array<double, 3> com = {0.0, 0.0, 0.0};
        double totalMass = 0.0;

        for (const auto& atom : atoms) {
            if (!atom) continue;
            double mass = atom->getMass();
            const auto& coor = atom->getCoor();
            for (int i = 0; i < 3; ++i) {
                com[i] += mass * coor[i];
            }
            totalMass += mass;
        }

        if (totalMass > 0.0) {
            for (double& x : com) x /= totalMass;
        }

        return com;
    }

    // Selection methods for CHARMM compatibility
    bool hasAtomType(const std::string& type) const {
        return std::any_of(atoms.begin(), atoms.end(),
            [&type](const std::shared_ptr<Atom>& atom) {
                return atom && atom->getType() == type;
            });
    }

    std::vector<std::shared_ptr<Atom>> selectAtoms(
        const std::function<bool(const Atom&)>& predicate) const {
        std::vector<std::shared_ptr<Atom>> selected;
        for (const auto& atom : atoms) {
            if (atom && predicate(*atom)) {
                selected.push_back(atom);
            }
        }
        return selected;
    }

    // Additional getters
    SecondaryStructure getSecondaryStructure() const noexcept { return secStruct; }
    const SheetStrand& getSheetInfo() const noexcept { return sheetInfo; }
    const SSBond& getSSBond() const noexcept { return ssbond; }

    // Additional setters
    void setSecondaryStructure(SecondaryStructure ss) { secStruct = ss; }
    void setSheetInfo(const SheetStrand& si) { sheetInfo = si; }
    void setSSBond(const SSBond& sb) { ssbond = sb; }

    // PDB format utilities
    std::string getResidueID() const {
        // Combine residue number and insertion code (e.g., "153A")
        if (inscode == ' ') {
            return std::to_string(ires);
        }
        return std::to_string(ires) + inscode;
    }

    void setResidueID(const std::string& resid) {
        // Parse residue ID (e.g., "153A" -> ires=153, inscode='A')
        size_t numLen = 0;
        try {
            ires = std::stoi(resid, &numLen);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid residue ID format");
        }
        
        if (numLen < resid.length()) {
            inscode = resid[numLen];
        } else {
            inscode = ' ';
        }
    }

    // Enhanced atom lookup methods
    std::shared_ptr<Atom> findAtomByPDBName(const std::string& pdbName) const {
        // Find atom by PDB formatted name
        return std::find_if(atoms.begin(), atoms.end(),
            [&pdbName](const std::shared_ptr<Atom>& atom) {
                return atom && atom->getFormattedAtomName() == pdbName;
            }) != atoms.end() ? *std::find_if(atoms.begin(), atoms.end(),
            [&pdbName](const std::shared_ptr<Atom>& atom) {
                return atom && atom->getFormattedAtomName() == pdbName;
            }) : nullptr;
    }

    void updateAtomMap() {
        atomMap.clear();
        for (const auto& atom : atoms) {
            if (atom) {
                // Store both raw and PDB-formatted names
                atomMap[atom->getType()] = atom;
                atomMap[atom->getFormattedAtomName()] = atom;
            }
        }
    }

    // CHARMM-style atom range
    std::pair<size_t, size_t> getAtomRange() const {
        return {0, atoms.size()};  // Equivalent to IBASE(IRES) to IBASE(IRES+1)
    }

private:
    // CHARMM standard fields
    std::string resname;     ///< Residue name (RESNAME)
    int ires;                ///< Residue number (IRES)
    std::string segid;       ///< Segment ID (SEGID)
    int iseg;                ///< Segment number (ISEG)
    int igro;                ///< Group number (IGRO)
    char chain;              ///< Chain identifier
    char inscode;            ///< Insertion code
    int move;                ///< Movement flag (MOVE)
    int ignore;              ///< Ignore flag (IGNORE)
    int constrain;           ///< Constraint flag (CONSTRAIN)

    // Atom storage
    std::vector<std::shared_ptr<Atom>> atoms;  ///< Atoms in residue
    std::unordered_map<std::string, std::shared_ptr<Atom>> atomMap;  ///< Quick atom lookup by type

    // Scalar properties (SCA1-9)
    std::array<double, 9> scalar;  ///< User-defined scalar properties

    // Secondary structure information
    SecondaryStructure secStruct{SecondaryStructure::NONE};
    SheetStrand sheetInfo;
    SSBond ssbond;
};

} // namespace data
} // namespace pygcmc

#endif // PYGCMC_DATA_RESIDUE_HPP