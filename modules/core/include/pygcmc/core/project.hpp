// modules/core/include/pygcmc/core/project.hpp

#pragma once

#include <string>
#include <memory>
#include <vector>
#include "pygcmc/core/system.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/ff_parser.hpp"

namespace pygcmc {
namespace core {

class Structure;
class ForceField;

class Project {
public:
    explicit Project(const std::string& name = "");
    ~Project();

    // Load structure from PDB and topology files
    std::shared_ptr<Structure> load_structure(const std::string& pdb_file, 
                                            const std::string& top_file);

    // Load force field from parameter files
    std::shared_ptr<ForceField> load_forcefield(const std::vector<std::string>& param_files);

    // Getters
    const std::string& get_name() const { return name_; }
    
private:
    std::string name_;
    std::vector<std::shared_ptr<Structure>> structures_;
    std::vector<std::shared_ptr<ForceField>> forcefields_;
}; 

} // namespace core
} // namespace pygcmc 