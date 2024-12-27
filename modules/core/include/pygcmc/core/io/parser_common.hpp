// modules/core/include/pygcmc/core/io/parser_common.hpp

#ifndef PYGCMC_CORE_IO_PARSER_COMMON_HPP
#define PYGCMC_CORE_IO_PARSER_COMMON_HPP

#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <array>
#include <map>
#include <tuple>
#include <unordered_map>
#include <functional>

// 现有的内容...

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief Base parser exception class
 */
class ParserError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief File format error
 */
class FormatError : public ParserError {
public:
    using ParserError::ParserError;
};

/**
 * @brief File access error
 */
class FileError : public ParserError {
public:
    using ParserError::ParserError;
};

/**
 * @brief Common atom structure for PDB files
 */
struct PDBAtom {
    int serial;              ///< Atom serial number
    std::string name;        ///< Atom name
    std::string residue;     ///< Residue name
    int sequence;            ///< Residue sequence number
    double x, y, z;          ///< Atomic coordinates
    double charge;           ///< Atomic charge
    std::string type;        ///< Atom type

    PDBAtom(int serial_ = 0, const std::string& name_ = "", 
            const std::string& residue_ = "", int sequence_ = 0,
            double x_ = 0.0, double y_ = 0.0, double z_ = 0.0,
            double charge_ = 0.0, const std::string& type_ = "")
        : serial(serial_), name(name_), residue(residue_),
          sequence(sequence_), x(x_), y(y_), z(z_),
          charge(charge_), type(type_) {}

    bool is_valid() const {
        return serial > 0 && !name.empty() && !residue.empty() &&
               std::isfinite(x) && std::isfinite(y) && std::isfinite(z) &&
               std::isfinite(charge) && !type.empty();
    }

    std::array<double, 3> position() const {
        return {x, y, z};
    }
};

/**
 * @brief ITP Atom structure
 */
struct ITPAtom {
    std::string name;    ///< Atom name
    std::string type;    ///< Atom type
    int resid;           ///< Residue ID
    std::string resname; ///< Residue name
    double charge;       ///< Atomic charge

    ITPAtom(const std::string& name_ = "", const std::string& type_ = "",
            int resid_ = 0, const std::string& resname_ = "", double charge_ = 0.0)
        : name(name_), type(type_), resid(resid_), resname(resname_), charge(charge_) {}

    bool is_valid() const {
        return !name.empty() && !type.empty() && resid > 0 && !resname.empty() && std::isfinite(charge);
    }
};

/**
 * @brief PSF Bond structure
 */
struct PSFBond {
    int atom1; ///< First atom index
    int atom2; ///< Second atom index

    PSFBond(int a1 = 0, int a2 = 0) : atom1(a1), atom2(a2) {}

    bool is_valid() const {
        return atom1 > 0 && atom2 > 0 && atom1 != atom2;
    }
};

/**
 * @brief PSF Topology structure
 */
struct PSFTopology {
    std::vector<PSFBond> bonds;

    bool is_valid() const {
        // 添加具体的验证逻辑，例如检查是否有重复的键等
        return true;
    }
};

/**
 * @brief Topology Atom Type structure
 */
struct TopAtomType {
    std::string name; ///< Atom type name
    int type;         ///< Atom type number
    double charge;    ///< Atomic charge
    double mass;      ///< Atomic mass

    TopAtomType(const std::string& name_ = "", int type_ = 0,
               double charge_ = 0.0, double mass_ = 0.0)
        : name(name_), type(type_), charge(charge_), mass(mass_) {}

    bool is_valid() const {
        return !name.empty() && type > 0 && std::isfinite(charge) && std::isfinite(mass);
    }
};

/**
 * @brief Topology structure
 */
struct Topology {
    std::vector<TopAtomType> atom_types;

    bool is_valid() const {
        // 添加具体的验证逻辑，例如检查是否有重复的原子类型等
        return true;
    }
};

/**
 * @brief Force field parameter pair structure
 */
struct ForceFieldPair {
    double param1;
    double param2;

    ForceFieldPair(double p1 = 0.0, double p2 = 0.0) 
        : param1(p1), param2(p2) {}
};

struct PairStringHash {
    std::size_t operator()(const std::pair<std::string, std::string>& p) const {
        return std::hash<std::string>()(p.first) ^ (std::hash<std::string>()(p.second) << 1);
    }
};

using NBMap = std::unordered_map<std::string, ForceFieldPair>;
using NBFixMap = std::unordered_map<std::pair<std::string, std::string>, ForceFieldPair, PairStringHash>;

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PARSER_COMMON_HPP
