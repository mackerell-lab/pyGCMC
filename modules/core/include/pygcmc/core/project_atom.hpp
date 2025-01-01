#pragma once

#include <memory>
#include "pygcmc/core/io/pdb_parser.hpp"

namespace pygcmc {
namespace core {

class ProjectAtom {
public:
    ProjectAtom();
    explicit ProjectAtom(const io::PDBAtom& atom);
    
    // Getters
    int get_serial() const;
    std::string get_name() const;
    std::string get_residue() const;
    int get_sequence() const;
    char get_chain() const;
    char get_alt_loc() const;
    char get_insertion_code() const;
    double get_x() const;
    double get_y() const;
    double get_z() const;
    double get_occupancy() const;
    double get_temp_factor() const;
    std::string get_element() const;
    std::string get_charge() const;
    std::string get_type() const;
    std::string get_topo_type() const;
    double get_topo_charge() const;
    double get_topo_mass() const;
    double get_forcefield_epsilon() const;
    double get_forcefield_rmin() const;

    // Setters
    void set_serial(int val);
    void set_name(const std::string& val);
    void set_residue(const std::string& val);
    void set_sequence(int val);
    void set_chain(char val);
    void set_alt_loc(char val);
    void set_insertion_code(char val);
    void set_x(double val);
    void set_y(double val);
    void set_z(double val);
    void set_occupancy(double val);
    void set_temp_factor(double val);
    void set_element(const std::string& val);
    void set_charge(const std::string& val);
    void set_type(const std::string& val);
    void set_topo_type(const std::string& val);
    void set_topo_charge(double val);
    void set_topo_mass(double val);
    void set_forcefield_epsilon(double val);
    void set_forcefield_rmin(double val);

    // Status checks
    bool is_valid() const;
    bool has_topology_info() const;
    bool has_forcefield_info() const;

    // Get underlying pointer
    std::shared_ptr<io::PDBAtom> get_ptr() const { return ptr; }

private:
    std::shared_ptr<io::PDBAtom> ptr;
};

} // namespace core
} // namespace pygcmc 