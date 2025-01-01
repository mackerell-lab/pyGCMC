// modules/core/src/project.cpp

#include "pygcmc/core/project.hpp"
#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace pygcmc {
namespace core {

Project::Project(const std::string& name) : name_(name) {}

Project::~Project() {
    // Clear vectors in reverse order of dependency
    forcefields_.clear();  // Clear force fields first
    structures_.clear();   // Then clear structures
}

Structure Project::load_structure(const std::string& pdb_file, const std::string& top_file) {
    Structure structure;
    
    // Parse PDB file
    auto [box_vec, residues] = io::PDBParser::parse(pdb_file);
    
    // Convert vector box to array box if present
    if (box_vec && box_vec->size() == 6) {
        std::array<double, 6> box_arr;
        std::copy(box_vec->begin(), box_vec->end(), box_arr.begin());
        structure.set_box(std::make_optional(box_arr));
    } else {
        structure.set_box(std::nullopt);
    }
    
    // Add residues to structure
    for (const auto& residue : residues) {
        // Create a new shared_ptr to a copy of the residue
        auto residue_ptr = std::make_shared<io::IOResidue>();
        *residue_ptr = residue;  // Use copy assignment

        // Create shared_ptr for each atom and update atom_ptrs
        residue_ptr->atom_ptrs.clear();  // Clear existing pointers
        for (const auto& atom : residue.atoms) {
            auto atom_ptr = std::make_shared<io::PDBAtom>(atom);
            residue_ptr->atom_ptrs.push_back(atom_ptr);
            structure.add_atom(atom_ptr);  // Add atom to structure's atoms_ vector
        }
        
        structure.add_residue(residue_ptr);
    }
    
    // If topology file is provided, load it
    if (!top_file.empty()) {
        // Create TopParser instance and parse the file with includes
        io::TopParser top_parser;
        if (!top_parser.parse_with_includes(top_file)) {
            throw std::runtime_error("Failed to parse topology file: " + top_file);
        }
        
        // Update all atoms with topology information using update_pdb_atoms
        std::vector<io::PDBAtom*> atom_ptrs;
        for (const auto& residue : structure.residues()) {
            for (const auto& atom : residue->atom_ptrs) {
                if (atom) {
                    atom_ptrs.push_back(atom.get());
                }
            }
        }
        
        int updated = top_parser.update_pdb_atoms(atom_ptrs);
        if (updated == 0) {
            throw std::runtime_error("No atoms were updated with topology information");
        }
    }
    
    // Store structure in project
    structures_.push_back(std::make_shared<Structure>(structure));
    structure_ = structures_.back();
    
    return structure;
}

std::shared_ptr<ForceField> Project::load_forcefield(const std::vector<std::string>& param_files) {
    // Create new force field
    auto ff = std::make_shared<ForceField>();
    
    try {
        // Parse each parameter file
        for (const auto& file : param_files) {
            io::FFParser parser;
            if (!parser.parse(file)) {
                throw std::runtime_error("Failed to parse force field file: " + file);
            }
            
            // Get parameters
            const auto& nonbonded = parser.get_nonbonded_params();
            const auto& nbfix = parser.get_nbfix_params();
            
            // Copy parameters into force field
            ff->nonbonded_params().insert(nonbonded.begin(), nonbonded.end());
            ff->nbfix_params().insert(nbfix.begin(), nbfix.end());

            // If we have a structure loaded, update its atoms with force field parameters
            if (structure_) {
                std::vector<io::PDBAtom*> atom_ptrs;
                for (const auto& residue : structure_->residues()) {
                    for (const auto& atom : residue->atom_ptrs) {
                        if (atom) {
                            atom_ptrs.push_back(atom.get());
                        }
                    }
                }
                parser.update_pdb_atoms(atom_ptrs);
            }
        }
        
        // Store force field in the project
        forcefields_.push_back(ff);
        forcefield_ = ff;  // Store current force field
        
        return ff;
        
    } catch (const std::exception& e) {
        // Clean up in case of error
        if (std::find(forcefields_.begin(), forcefields_.end(), ff) != forcefields_.end()) {
            forcefields_.erase(
                std::remove(forcefields_.begin(), forcefields_.end(), ff),
                forcefields_.end()
            );
        }
        forcefield_ = nullptr;  // Clear current force field
        throw;
    }
}

void Project::print_atom_info(const ProjectAtom& atom) const {
    std::cout << "Debug info for atom " << atom.get_name() << " in residue " 
              << atom.get_residue() << " " << atom.get_sequence() << ":\n";
    std::cout << "  Basic info:\n";
    std::cout << "    serial: " << atom.get_serial() << "\n";
    std::cout << "    name: " << atom.get_name() << "\n";
    std::cout << "    residue: " << atom.get_residue() << "\n";
    std::cout << "    sequence: " << atom.get_sequence() << "\n";
    std::cout << "    position: (" << atom.get_x() << ", " << atom.get_y() << ", " << atom.get_z() << ")\n";
    std::cout << "  Topology info:\n";
    std::cout << "    type: " << atom.get_type() << "\n";
    std::cout << "    topo_type: " << atom.get_topo_type() << "\n";
    std::cout << "    topo_charge: " << atom.get_topo_charge() << "\n";
    std::cout << "    topo_mass: " << atom.get_topo_mass() << "\n";
    std::cout << "  Force field info:\n";
    std::cout << "    forcefield_epsilon: " << atom.get_forcefield_epsilon() << "\n";
    std::cout << "    forcefield_rmin: " << atom.get_forcefield_rmin() << "\n";
    std::cout << "  Status:\n";
    std::cout << "    has_topology_info: " << (atom.has_topology_info() ? "True" : "False") << "\n";
    std::cout << "    has_forcefield_info: " << (atom.has_forcefield_info() ? "True" : "False") << "\n\n";
}

void Project::print_detailed_atom_info(const ProjectAtom& atom) const {
    std::cout << std::left << std::setw(10) << atom.get_residue()
              << std::setw(6) << atom.get_sequence()
              << std::setw(5) << atom.get_name()
              << std::setw(8) << atom.get_type()
              << std::setw(10) << atom.get_topo_type()
              << std::fixed << std::setprecision(3)
              << std::setw(10) << atom.get_topo_charge()
              << std::setw(10) << atom.get_topo_mass()
              << std::setw(10) << atom.get_forcefield_epsilon()
              << std::setw(10) << atom.get_forcefield_rmin()
              << "\n";
}

void Project::print_atom_table_header() const {
    std::cout << "\nDetailed Atom Information:\n";
    std::cout << std::string(120, '-') << "\n";
    std::cout << std::left 
              << std::setw(10) << "Residue"
              << std::setw(6) << "Seq"
              << std::setw(5) << "Name"
              << std::setw(8) << "Type"
              << std::setw(10) << "TopoType"
              << std::setw(10) << "Charge"
              << std::setw(10) << "Mass"
              << std::setw(10) << "Epsilon"
              << std::setw(10) << "Rmin"
              << "\n";
    std::cout << std::string(120, '-') << "\n";
}

void Project::print_all_atoms() const {
    if (!structure_) {
        std::cout << "No structure loaded.\n";
        return;
    }

    print_atom_table_header();
    for (const auto& atom_ptr : structure_->atoms()) {
        if (atom_ptr) {
            ProjectAtom proj_atom(*atom_ptr);  // Convert PDBAtom to ProjectAtom
            print_detailed_atom_info(proj_atom);
        }
    }
    std::cout << "\n";
}

void Project::print_forcefield_info() const {
    if (!forcefield_) {
        std::cout << "No force field loaded." << std::endl;
        return;
    }

    std::cout << "\nNonbonded Parameters Loaded:" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(10) << "Type" 
              << std::setw(15) << "Epsilon" 
              << std::setw(15) << "Rmin" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    for (const auto& [type, params] : forcefield_->nonbonded_params()) {
        std::cout << std::left << std::setw(10) << type 
                 << std::setw(15) << params.epsilon 
                 << std::setw(15) << params.rmin << std::endl;
    }
    std::cout << std::endl;

    std::cout << "\n=== Force Field Parameters Summary ===\n\n";
    print_global_parameters();
    print_nbfix_info();
}

void Project::print_global_parameters() const {
    if (!forcefield_) return;

    std::cout << "Global Nonbonded Parameters:\n";
    std::cout << "  cutoff: " << forcefield_->get_cutoff() << " Å\n";
    std::cout << "  switching: " << forcefield_->get_switching() << " Å\n";
    std::cout << "  pairlist_distance: " << forcefield_->get_pairlist_distance() << " Å\n\n";
}

void Project::print_nbfix_info() const {
    if (!forcefield_) {
        std::cout << "No force field loaded." << std::endl;
        return;
    }

    std::cout << "\nNBFIX Parameters:" << std::endl;
    std::cout << "Total NBFIX parameters: " << forcefield_->nbfix_params().size() << std::endl;
    std::cout << "\nAll NBFIX parameters:" << std::endl;

    for (const auto& [types, params] : forcefield_->nbfix_params()) {
        std::cout << "  " << types.first << "-" << types.second 
                 << ": epsilon=" << params.epsilon 
                 << ", rmin=" << params.rmin << std::endl;
    }
    std::cout << std::endl;
}

} // namespace core
} // namespace pygcmc