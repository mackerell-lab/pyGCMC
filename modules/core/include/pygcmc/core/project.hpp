// modules/core/include/pygcmc/core/project.hpp

#pragma once

#include <string>
#include <memory>
#include <vector>
#include <optional>
#include "pygcmc/core/system.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/ff_parser.hpp"
#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"
#include "pygcmc/core/project_atom.hpp"
#include "pygcmc/core/project_residue.hpp"

namespace pygcmc {
namespace core {

class Project {
public:
    explicit Project(const std::string& name);
    ~Project();

    // Get project name
    std::string get_name() const { return name_; }

    // Load structure from PDB and optionally TOP files
    Structure load_structure(const std::string& pdb_file, const std::string& top_file = "");

    // Create a new empty structure
    Structure create_structure();

    // Load force field from parameter files
    std::shared_ptr<ForceField> load_forcefield(const std::vector<std::string>& param_files);

    // Methods for detailed information
    void print_atom_info(const ProjectAtom& atom) const;
    void print_detailed_atom_info(const ProjectAtom& atom) const;
    void print_atom_table_header() const;
    void print_all_atoms() const;
    void print_forcefield_info() const;
    void print_nbfix_info() const;
    void print_global_parameters() const;

private:
    std::string name_;
    std::shared_ptr<Structure> structure_;
    std::shared_ptr<ForceField> forcefield_;
    std::vector<std::shared_ptr<Structure>> structures_;
    std::vector<std::shared_ptr<ForceField>> forcefields_;
}; 

} // namespace core
} // namespace pygcmc 