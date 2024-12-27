// modules/core/include/pygcmc/core/io/parser.hpp

#ifndef PYGCMC_CORE_IO_PARSER_HPP
#define PYGCMC_CORE_IO_PARSER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include "pygcmc/core/force_field.hpp"

namespace pygcmc {
namespace core {
namespace io {

// Add base parser exception
class ParserError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

// ---------------------------
// 1. PDB Parser
// ---------------------------

struct PDBAtom {
    int serial;
    std::string name;
    std::string residue;
    int sequence;
    double x, y, z;
    double charge;
    int type;

    PDBAtom(int serial_ = 0, const std::string& name_ = "", 
            const std::string& residue_ = "", int sequence_ = 0,
            double x_ = 0.0, double y_ = 0.0, double z_ = 0.0,
            double charge_ = 0.0, int type_ = 0)
        : serial(serial_), name(name_), residue(residue_),
          sequence(sequence_), x(x_), y(y_), z(z_),
          charge(charge_), type(type_) {}

    bool is_valid() const {
        return serial > 0 && !name.empty() && !residue.empty() &&
               std::isfinite(x) && std::isfinite(y) && std::isfinite(z) &&
               std::isfinite(charge);
    }
};

class PDBParser {
public:
    static std::pair<std::vector<double>, std::vector<PDBAtom>> parse(const std::string& filename);

protected:
    static bool parse_cryst1_line(const std::string& line, std::vector<double>& cell_params);
    static bool parse_atom_line(const std::string& line, PDBAtom& atom);
    static bool validate_pdb_structure(const std::vector<PDBAtom>& atoms);
};

// ---------------------------
// 2. PSF Parser
// ---------------------------

struct PSFBond {
    int atom1;
    int atom2;

    PSFBond(int a1 = 0, int a2 = 0) : atom1(a1), atom2(a2) {}

    bool is_valid() const {
        return atom1 > 0 && atom2 > 0 && atom1 != atom2;
    }
};

struct PSFTopology {
    std::vector<PSFBond> bonds;

    bool is_valid() const {
        return std::all_of(bonds.begin(), bonds.end(), 
                          [](const PSFBond& b) { return b.is_valid(); });
    }
};

class PSFParser {
public:
    static PSFTopology parse(const std::string& filename);

protected:
    static bool parse_bonds_section(std::istream& is, PSFTopology& topology);
};

// ---------------------------
// 3. TOP Parser
// ---------------------------

struct TopAtomType {
    std::string name;
    int type;
    double charge;
    double mass;

    TopAtomType(const std::string& name_ = "", int type_ = 0,
                double charge_ = 0.0, double mass_ = 0.0)
        : name(name_), type(type_), charge(charge_), mass(mass_) {}

    bool is_valid() const {
        return !name.empty() && type >= 0 && 
               std::isfinite(charge) && mass > 0.0;
    }
};

struct Topology {
    std::vector<TopAtomType> atom_types;

    bool is_valid() const {
        return std::all_of(atom_types.begin(), atom_types.end(),
                          [](const TopAtomType& at) { return at.is_valid(); });
    }
};

class TopParser {
public:
    static Topology parse(const std::string& filename);

protected:
    static bool parse_atomtypes_section(std::istream& is, Topology& top);
};

// ---------------------------
// 4. ITP Parser
// ---------------------------

struct ITPAtom {
    std::string name;
    std::string type;
    int resid;
    std::string resname;
    double charge;

    ITPAtom(const std::string& name_ = "", const std::string& type_ = "",
            int resid_ = 0, const std::string& resname_ = "", double charge_ = 0.0)
        : name(name_), type(type_), resid(resid_),
          resname(resname_), charge(charge_) {}

    bool is_valid() const {
        return !name.empty() && !type.empty() && 
               resid > 0 && !resname.empty() && 
               std::isfinite(charge);
    }
};

class ITPParser {
public:
    static std::vector<ITPAtom> parse(const std::string& filename);

protected:
    static bool parse_atoms_section(std::istream& is, std::vector<ITPAtom>& atoms);
};

// ---------------------------
// 5. Force Field Parser
// ---------------------------

class FFParser {
public:
    using NBMap = core::NBMap;
    using NBFixMap = core::NBFixMap;
    
    static std::pair<NBMap, NBFixMap> parse(const std::string& filename);

protected:
    static bool parse_nb_section(std::istream& is, NBMap& nb_params);
    static bool parse_nbfix_section(std::istream& is, NBFixMap& nbfix_params);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PARSER_HPP
