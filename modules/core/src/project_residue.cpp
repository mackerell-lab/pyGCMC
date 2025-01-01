#include "pygcmc/core/project_residue.hpp"

namespace pygcmc {
namespace core {

ProjectResidue::ProjectResidue() : ptr(std::make_shared<io::IOResidue>()) {}

ProjectResidue::ProjectResidue(const std::string& name, int seq_num, char chain)
    : ptr(std::make_shared<io::IOResidue>(name, seq_num, chain)) {}

std::string ProjectResidue::get_name() const { return ptr->name; }
int ProjectResidue::get_sequence_number() const { return ptr->sequence_number; }
char ProjectResidue::get_chain_id() const { return ptr->chain_id; }

std::vector<ProjectAtom> ProjectResidue::get_atoms() const {
    std::vector<ProjectAtom> result;
    for (const auto& atom_ptr : ptr->atom_ptrs) {
        if (atom_ptr) {
            result.emplace_back(*atom_ptr);
        }
    }
    return result;
}

void ProjectResidue::set_name(const std::string& val) { ptr->name = val; }
void ProjectResidue::set_sequence_number(int val) { ptr->sequence_number = val; }
void ProjectResidue::set_chain_id(char val) { ptr->chain_id = val; }

void ProjectResidue::set_atoms(const std::vector<ProjectAtom>& atoms) {
    ptr->atom_ptrs.clear();
    for (const auto& atom : atoms) {
        ptr->atom_ptrs.push_back(atom.get_ptr());
    }
}

std::array<double, 3> ProjectResidue::center_of_mass() const {
    return ptr->center_of_mass();
}

size_t ProjectResidue::atom_count() const {
    return ptr->atom_count();
}

} // namespace core
} // namespace pygcmc 