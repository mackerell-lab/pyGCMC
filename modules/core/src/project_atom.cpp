#include "pygcmc/core/project_atom.hpp"
#include <cmath>

namespace pygcmc {
namespace core {

ProjectAtom::ProjectAtom() : ptr(std::make_shared<io::PDBAtom>()) {}

ProjectAtom::ProjectAtom(const io::PDBAtom& atom) : ptr(std::make_shared<io::PDBAtom>(atom)) {}

// Getters
int ProjectAtom::get_serial() const { return ptr->serial; }
std::string ProjectAtom::get_name() const { return ptr->name; }
std::string ProjectAtom::get_residue() const { return ptr->residue; }
int ProjectAtom::get_sequence() const { return ptr->sequence; }
char ProjectAtom::get_chain() const { return ptr->chain; }
char ProjectAtom::get_alt_loc() const { return ptr->alt_loc; }
char ProjectAtom::get_insertion_code() const { return ptr->insertion_code; }
double ProjectAtom::get_x() const { return ptr->x; }
double ProjectAtom::get_y() const { return ptr->y; }
double ProjectAtom::get_z() const { return ptr->z; }
double ProjectAtom::get_occupancy() const { return ptr->occupancy; }
double ProjectAtom::get_temp_factor() const { return ptr->temp_factor; }
std::string ProjectAtom::get_element() const { return ptr->element; }
std::string ProjectAtom::get_charge() const { return ptr->charge; }
std::string ProjectAtom::get_type() const { return ptr->type; }
std::string ProjectAtom::get_topo_type() const { return ptr->topo_type; }
double ProjectAtom::get_topo_charge() const { return ptr->topo_charge; }
double ProjectAtom::get_topo_mass() const { return ptr->topo_mass; }
double ProjectAtom::get_forcefield_epsilon() const { return ptr->forcefield_epsilon; }
double ProjectAtom::get_forcefield_rmin() const { return ptr->forcefield_rmin; }

// Setters
void ProjectAtom::set_serial(int val) { ptr->serial = val; }
void ProjectAtom::set_name(const std::string& val) { ptr->name = val; }
void ProjectAtom::set_residue(const std::string& val) { ptr->residue = val; }
void ProjectAtom::set_sequence(int val) { ptr->sequence = val; }
void ProjectAtom::set_chain(char val) { ptr->chain = val; }
void ProjectAtom::set_alt_loc(char val) { ptr->alt_loc = val; }
void ProjectAtom::set_insertion_code(char val) { ptr->insertion_code = val; }
void ProjectAtom::set_x(double val) { ptr->x = val; }
void ProjectAtom::set_y(double val) { ptr->y = val; }
void ProjectAtom::set_z(double val) { ptr->z = val; }
void ProjectAtom::set_occupancy(double val) { ptr->occupancy = val; }
void ProjectAtom::set_temp_factor(double val) { ptr->temp_factor = val; }
void ProjectAtom::set_element(const std::string& val) { ptr->element = val; }
void ProjectAtom::set_charge(const std::string& val) { ptr->charge = val; }
void ProjectAtom::set_type(const std::string& val) { ptr->type = val; }
void ProjectAtom::set_topo_type(const std::string& val) { ptr->topo_type = val; }
void ProjectAtom::set_topo_charge(double val) { ptr->topo_charge = val; }
void ProjectAtom::set_topo_mass(double val) { ptr->topo_mass = val; }
void ProjectAtom::set_forcefield_epsilon(double val) { ptr->forcefield_epsilon = val; }
void ProjectAtom::set_forcefield_rmin(double val) { ptr->forcefield_rmin = val; }

// Status checks
bool ProjectAtom::is_valid() const { return ptr->is_valid(); }

bool ProjectAtom::has_topology_info() const {
    return !ptr->topo_type.empty() && !std::isnan(ptr->topo_charge) && !std::isnan(ptr->topo_mass);
}

bool ProjectAtom::has_forcefield_info() const {
    return !std::isnan(ptr->forcefield_epsilon) && !std::isnan(ptr->forcefield_rmin);
}

} // namespace core
} // namespace pygcmc 