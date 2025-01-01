#pragma once

#include <memory>
#include <vector>
#include <array>
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/project_atom.hpp"

namespace pygcmc {
namespace core {

class ProjectResidue {
public:
    ProjectResidue();
    ProjectResidue(const std::string& name, int seq_num = 0, char chain = ' ');

    // Getters
    std::string get_name() const;
    int get_sequence_number() const;
    char get_chain_id() const;
    std::vector<ProjectAtom> get_atoms() const;

    // Setters
    void set_name(const std::string& val);
    void set_sequence_number(int val);
    void set_chain_id(char val);
    void set_atoms(const std::vector<ProjectAtom>& atoms);

    // Methods
    std::array<double, 3> center_of_mass() const;
    size_t atom_count() const;

    // Get underlying pointer
    std::shared_ptr<io::IOResidue> get_ptr() const { return ptr; }

private:
    std::shared_ptr<io::IOResidue> ptr;
};

} // namespace core
} // namespace pygcmc 