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
    auto [box, residues] = io::PDBParser::parse(pdb_file);
    structure.set_box(box);
    
    // Add residues to structure
    for (const auto& residue : residues) {
        // Create a new shared_ptr to a copy of the residue
        auto residue_ptr = std::make_shared<io::IOResidue>();
        *residue_ptr = residue;  // Use copy assignment
        structure.add_residue(residue_ptr);
    }
    
    // If topology file is provided, load it
    if (!top_file.empty()) {
        // Create TopParser instance and parse the file
        io::TopParser top_parser;
        if (!top_parser.parse(top_file)) {
            throw std::runtime_error("Failed to parse topology file: " + top_file);
        }
        
        // Apply topology to each residue's atoms
        for (const auto& residue : structure.residues()) {
            for (const auto& atom : residue->atom_ptrs) {
                if (atom) {
                    // Get topology information for this atom
                    double charge, mass;
                    if (top_parser.get_atom_properties(atom->residue, atom->name, charge, mass)) {
                        atom->topo_charge = charge;
                        atom->topo_mass = mass;
                    }
                }
            }
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
            parser.parse(file);
            
            // Get parameters
            const auto& nonbonded = parser.get_nonbonded_params();
            const auto& nbfix = parser.get_nbfix_params();
            
            // Copy parameters into force field
            ff->nonbonded_params().insert(nonbonded.begin(), nonbonded.end());
            ff->nbfix_params().insert(nbfix.begin(), nbfix.end());
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
        std::cout << "No force field loaded.\n";
        return;
    }

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
    if (!forcefield_) return;

    const auto& nbfix = forcefield_->nbfix_params();
    std::cout << "NBFIX Parameters:\n";
    std::cout << "Total NBFIX parameters: " << nbfix.size() << "\n\n";

    if (!nbfix.empty()) {
        std::cout << "All NBFIX parameters:\n";
        for (const auto& [types, params] : nbfix) {
            std::cout << "  " << types.first << "-" << types.second 
                     << ": epsilon=" << params.epsilon 
                     << ", rmin=" << params.rmin << "\n";
        }
        std::cout << "\n";
    }
}

} // namespace core
} // namespace pygcmc