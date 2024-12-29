// modules/core/include/pygcmc/core/io/pdb_atom.hpp

#ifndef PYGCMC_CORE_IO_PDB_ATOM_HPP
#define PYGCMC_CORE_IO_PDB_ATOM_HPP

#include <string>
#include <cmath>

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief Represents an atom parsed from a PDB file.
 */
struct PDBAtom {
    int serial;               ///< Atom serial number
    std::string name;         ///< Atom name
    std::string residue;      ///< Residue name
    int sequence;             ///< Residue sequence number
    char chain;               ///< Chain identifier
    char alt_loc;             ///< Alternate location indicator
    char insertion_code;      ///< Insertion code
    double x, y, z;           ///< Atomic coordinates
    double occupancy;         ///< Occupancy
    double temp_factor;       ///< Temperature factor
    std::string element;      ///< Element symbol
    std::string charge;       ///< Charge
    std::string type;         ///< Atom type

    // Constructors
    PDBAtom()
        : serial(0), sequence(0), chain(' '), alt_loc(' '), insertion_code(' '),
          x(0.0), y(0.0), z(0.0), occupancy(0.0), temp_factor(0.0) {}

    PDBAtom(int serial_, const std::string& name_, const std::string& residue_,
            int sequence_, char chain_, char alt_loc_, char insertion_code_,
            double x_, double y_, double z_, double occupancy_,
            double temp_factor_, const std::string& element_,
            const std::string& charge_, const std::string& type_)
        : serial(serial_), name(name_), residue(residue_), sequence(sequence_),
          chain(chain_), alt_loc(alt_loc_), insertion_code(insertion_code_),
          x(x_), y(y_), z(z_), occupancy(occupancy_),
          temp_factor(temp_factor_), element(element_), charge(charge_), type(type_) {}

    bool is_valid() const {
        // 基本要求：有效的序号、名称、残基名和序列号
        bool basic_valid = serial > 0 && !name.empty() && !residue.empty() && sequence > 0;
        
        // 坐标必须是有限数
        bool coords_valid = std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
        
        // 元素和类型至少有一个不为空
        bool type_valid = !element.empty() || !type.empty();
        
        return basic_valid && coords_valid && type_valid;
    }
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PDB_ATOM_HPP
