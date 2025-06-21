#pragma once

#ifndef PYGCMC_MODEL_RESIDUE_COMPOSITE_HPP
#define PYGCMC_MODEL_RESIDUE_COMPOSITE_HPP

#include "../common/ModelInterface.hpp"
#include "../common/ModelConstants.hpp"
#include "../atom/AtomMain.hpp"
#include <vector>
#include <memory>
#include <string>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <functional>

namespace pygcmc {
namespace model {
namespace residue {

/**
 * @brief Core Residue class that handles atom composition
 * This class manages a collection of atoms that form a residue
 */
class ResidueComposite : public common::IValidatable, public common::IIdentifiable {
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
    ResidueComposite() : 
        resname(""), ires(0), segid(""), iseg(0), igro(0),
        chain(common::constants::DEFAULT_CHAIN),
        inscode(common::constants::DEFAULT_INSCODE),
        move(common::constants::DEFAULT_MOVE),
        ignore(common::constants::DEFAULT_IGNORE),
        constrain(common::constants::DEFAULT_CONSTRAIN),
        hetatm(false), com{0.0, 0.0, 0.0}
    {}

    ResidueComposite(const std::string& resname, int ires, 
                    const std::string& segid = "", int iseg = 0,
                    char chain = common::constants::DEFAULT_CHAIN, 
                    char inscode = common::constants::DEFAULT_INSCODE) :
        resname(resname), ires(ires), segid(segid), iseg(iseg), igro(0),
        chain(chain), inscode(inscode),
        move(common::constants::DEFAULT_MOVE),
        ignore(common::constants::DEFAULT_IGNORE),
        constrain(common::constants::DEFAULT_CONSTRAIN),
        hetatm(false), com{0.0, 0.0, 0.0}
    {}

    virtual ~ResidueComposite() = default;

    // IIdentifiable interface
    int get_id() const override { return ires; }
    void set_id(int id) override { 
        if (id <= 0) throw std::invalid_argument("Invalid residue number");
        ires = id; 
    }

    // IValidatable interface
    bool is_valid() const override {
        return !resname.empty() && ires > 0 && !segid.empty() &&
               std::all_of(atoms.begin(), atoms.end(),
                          [](const std::shared_ptr<atom::Atom>& atom) {
                              return atom && atom->is_valid();
                          });
    }

    // Core getters
    const std::string& get_resname() const noexcept { return resname; }
    int get_ires() const noexcept { return ires; }
    const std::string& get_segid() const noexcept { return segid; }
    int get_iseg() const noexcept { return iseg; }
    char get_chain() const noexcept { return chain; }
    char get_inscode() const noexcept { return inscode; }
    bool is_hetatm() const noexcept { return hetatm; }

    // Core setters
    void set_resname(const std::string& name) { 
        if (name.empty()) throw std::invalid_argument("Empty residue name");
        resname = name; 
    }
    void set_chain(char ch) { chain = ch; }
    void set_inscode(char ins) { inscode = ins; }
    void set_hetatm(bool het) { hetatm = het; }

    // Atom management - core composition logic
    void add_atom(const atom::Atom& atom) {
        // Verify atom belongs to this residue
        if (atom.get_resname() != resname || atom.get_ires() != ires ||
            atom.get_segid() != segid || atom.get_iseg() != iseg) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        atoms.push_back(std::make_shared<atom::Atom>(atom));
        update_atom_map();
    }

    void add_atom(std::shared_ptr<atom::Atom> atom) {
        if (!atom) throw std::invalid_argument("Null atom pointer");
        if (atom->get_resname() != resname || atom->get_ires() != ires ||
            atom->get_segid() != segid || atom->get_iseg() != iseg) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        atoms.push_back(atom);
        update_atom_map();
    }

    void remove_atom(const std::string& type) {
        atoms.erase(std::remove_if(atoms.begin(), atoms.end(),
            [&type](const std::shared_ptr<atom::Atom>& atom) {
                return atom && atom->get_type() == type;
            }), atoms.end());
        update_atom_map();
    }

    void clear_atoms() {
        atoms.clear();
        atomMap.clear();
    }

    // Atom access
    const std::vector<std::shared_ptr<atom::Atom>>& get_atoms() const noexcept { 
        return atoms; 
    }

    std::shared_ptr<atom::Atom> find_atom(const std::string& type) const {
        auto it = atomMap.find(type);
        return (it != atomMap.end()) ? it->second : nullptr;
    }

    size_t atom_count() const noexcept {
        return atoms.size();
    }

    // Atom selection
    std::vector<std::shared_ptr<atom::Atom>> select_atoms(
        const std::function<bool(const atom::Atom&)>& predicate) const {
        std::vector<std::shared_ptr<atom::Atom>> selected;
        for (const auto& atom : atoms) {
            if (atom && predicate(*atom)) {
                selected.push_back(atom);
            }
        }
        return selected;
    }

    bool has_atom_type(const std::string& type) const {
        return atomMap.find(type) != atomMap.end();
    }

    // Center of mass calculation
    void calculate_center_of_mass() {
        com = {0.0, 0.0, 0.0};
        double totalMass = 0.0;

        for (const auto& atom : atoms) {
            if (!atom) continue;
            double mass = atom->get_mass();
            const auto& coor = atom->get_coor();
            for (int i = 0; i < 3; ++i) {
                com[i] += mass * coor[i];
            }
            totalMass += mass;
        }

        if (totalMass > 0.0) {
            for (double& x : com) x /= totalMass;
        }
    }

    const std::array<double, 3>& get_center_of_mass() const noexcept {
        return com;
    }

    // Secondary structure
    SecondaryStructure get_secondary_structure() const noexcept { return secStruct; }
    void set_secondary_structure(SecondaryStructure ss) { secStruct = ss; }

    const SheetStrand& get_sheet_info() const noexcept { return sheetInfo; }
    void set_sheet_info(const SheetStrand& si) { sheetInfo = si; }

    const SSBond& get_ssbond() const noexcept { return ssbond; }
    void set_ssbond(const SSBond& sb) { ssbond = sb; }

    // Utility methods
    std::pair<size_t, size_t> get_atom_range() const {
        return {0, atoms.size()};
    }

protected:
    void update_atom_map() {
        atomMap.clear();
        for (const auto& atom : atoms) {
            if (atom) {
                atomMap[atom->get_type()] = atom;
            }
        }
    }

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
    bool hetatm;             ///< HETATM flag

    // Atom storage
    std::vector<std::shared_ptr<atom::Atom>> atoms;  ///< Atoms in residue
    std::unordered_map<std::string, std::shared_ptr<atom::Atom>> atomMap;  ///< Quick lookup

    // Secondary structure information
    SecondaryStructure secStruct{SecondaryStructure::NONE};
    SheetStrand sheetInfo;
    SSBond ssbond;

    // Scalar properties
    std::array<double, 9> scalar;  ///< User-defined scalar properties
    std::array<double, 3> com;     ///< Center of mass
};

} // namespace residue
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_RESIDUE_COMPOSITE_HPP 