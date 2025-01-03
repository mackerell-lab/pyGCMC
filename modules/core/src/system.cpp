// modules/core/src/system.cpp

#include "pygcmc/core/system.hpp"
#include "pygcmc/core/project.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <variant>

namespace pygcmc {
namespace core {

// Residue management
size_t System::add_residue(const std::string& name) {
    residues_.emplace_back(name);
    return residues_.size() - 1;
}

void System::remove_residue(size_t index) {
    validate_residue_index(index);
    residues_.erase(residues_.begin() + index);
}

size_t System::get_residue_count() const {
    return residues_.size();
}

const Residue& System::get_residue(size_t index) const {
    validate_residue_index(index);
    return residues_[index];
}

Residue& System::get_residue(size_t index) {
    validate_residue_index(index);
    return residues_[index];
}

// Particle management
size_t System::add_particle(size_t residue_index, const Particle& particle) {
    validate_residue_index(residue_index);
    residues_[residue_index].particles.push_back(particle);
    return residues_[residue_index].particles.size() - 1;
}

void System::remove_particle(size_t residue_index, size_t particle_index) {
    validate_particle_index(residue_index, particle_index);
    residues_[residue_index].particles.erase(
        residues_[residue_index].particles.begin() + particle_index);
}

size_t System::get_particle_count(size_t residue_index) const {
    validate_residue_index(residue_index);
    return residues_[residue_index].particles.size();
}

const Particle& System::get_particle(size_t residue_index, size_t particle_index) const {
    validate_particle_index(residue_index, particle_index);
    return residues_[residue_index].particles[particle_index];
}

Particle& System::get_particle(size_t residue_index, size_t particle_index) {
    validate_particle_index(residue_index, particle_index);
    return residues_[residue_index].particles[particle_index];
}

// Mass management
double System::get_particle_mass(size_t residue_index, size_t particle_index) const {
    validate_particle_index(residue_index, particle_index);
    return residues_[residue_index].particles[particle_index].mass;
}

void System::set_particle_mass(size_t residue_index, size_t particle_index, double mass) {
    validate_particle_index(residue_index, particle_index);
    if (mass < 0.0) {
        throw std::invalid_argument("Particle mass cannot be negative");
    }
    residues_[residue_index].particles[particle_index].mass = mass;
    if (mass == 0.0) {
        residues_[residue_index].particles[particle_index].is_virtual = true;
    }
}

// Virtual site management
void System::set_virtual_site(size_t residue_index, size_t particle_index, bool is_virtual) {
    validate_particle_index(residue_index, particle_index);
    residues_[residue_index].particles[particle_index].is_virtual = is_virtual;
    if (is_virtual) {
        residues_[residue_index].particles[particle_index].mass = 0.0;
    }
}

bool System::is_virtual_site(size_t residue_index, size_t particle_index) const {
    validate_particle_index(residue_index, particle_index);
    return residues_[residue_index].particles[particle_index].is_virtual;
}

// Constraint management
size_t System::add_constraint(size_t residue1, size_t particle1,
                            size_t residue2, size_t particle2, double distance) {
    validate_particle_index(residue1, particle1);
    validate_particle_index(residue2, particle2);
    if (distance <= 0.0) {
        throw std::invalid_argument("Constraint distance must be positive");
    }
    if (get_particle_mass(residue1, particle1) == 0.0 ||
        get_particle_mass(residue2, particle2) == 0.0) {
        throw std::invalid_argument("Cannot constrain virtual sites");
    }

    // Calculate global particle indices
    size_t global_particle1 = particle1;
    size_t global_particle2 = particle2;
    for (size_t i = 0; i < residue1; ++i) {
        global_particle1 += residues_[i].particles.size();
    }
    for (size_t i = 0; i < residue2; ++i) {
        global_particle2 += residues_[i].particles.size();
    }

    constraints_.emplace_back(global_particle1, global_particle2, distance);
    return constraints_.size() - 1;
}

void System::remove_constraint(size_t index) {
    validate_constraint_index(index);
    constraints_.erase(constraints_.begin() + index);
}

size_t System::get_constraint_count() const {
    return constraints_.size();
}

const Constraint& System::get_constraint(size_t index) const {
    validate_constraint_index(index);
    return constraints_[index];
}

void System::get_constraint_parameters(size_t index, size_t& residue1, size_t& particle1,
                                     size_t& residue2, size_t& particle2, double& distance) const {
    validate_constraint_index(index);
    const auto& constraint = constraints_[index];
    distance = constraint.distance;

    // Find residues and local particle indices
    size_t total_particles = 0;
    size_t global_particle1 = constraint.particle1;
    size_t global_particle2 = constraint.particle2;
    bool found1 = false, found2 = false;

    for (size_t i = 0; i < residues_.size(); ++i) {
        size_t residue_size = residues_[i].particles.size();
        
        if (!found1 && global_particle1 < total_particles + residue_size) {
            residue1 = i;
            particle1 = global_particle1 - total_particles;
            found1 = true;
        }
        
        if (!found2 && global_particle2 < total_particles + residue_size) {
            residue2 = i;
            particle2 = global_particle2 - total_particles;
            found2 = true;
        }
        
        total_particles += residue_size;
        if (found1 && found2) break;
    }

    if (!found1 || !found2) {
        throw std::runtime_error("Constraint references non-existent particles");
    }
}

void System::set_constraint_parameters(size_t index, size_t residue1, size_t particle1,
                                     size_t residue2, size_t particle2, double distance) {
    validate_constraint_index(index);
    validate_particle_index(residue1, particle1);
    validate_particle_index(residue2, particle2);
    if (distance <= 0.0) {
        throw std::invalid_argument("Constraint distance must be positive");
    }
    if (get_particle_mass(residue1, particle1) == 0.0 ||
        get_particle_mass(residue2, particle2) == 0.0) {
        throw std::invalid_argument("Cannot constrain virtual sites");
    }
    constraints_[index] = Constraint(particle1, particle2, distance);
}

// Force management
void System::add_force(std::shared_ptr<Force> force) {
    if (!force) {
        throw std::invalid_argument("Cannot add null force");
    }
    forces_.push_back(force);
}

void System::remove_force(size_t index) {
    validate_force_index(index);
    forces_.erase(forces_.begin() + index);
}

size_t System::get_force_count() const {
    return forces_.size();
}

std::shared_ptr<Force> System::get_force(size_t index) const {
    validate_force_index(index);
    return forces_[index];
}

// Periodic boundary conditions
void System::set_periodic_box_vectors(const std::array<double, 3>& a,
                                    const std::array<double, 3>& b,
                                    const std::array<double, 3>& c) {
    validate_box_vectors(a, b, c);
    box_vectors_[0] = a;
    box_vectors_[1] = b;
    box_vectors_[2] = c;
    has_periodic_boundary_ = true;
}

void System::get_periodic_box_vectors(std::array<double, 3>& a,
                                    std::array<double, 3>& b,
                                    std::array<double, 3>& c) const {
    a = box_vectors_[0];
    b = box_vectors_[1];
    c = box_vectors_[2];
}

bool System::uses_periodic_boundary_conditions() const {
    if (!has_periodic_boundary_) return false;
    for (const auto& force : forces_) {
        if (force->uses_periodic_boundary_conditions()) {
            return true;
        }
    }
    return false;
}

// System state
double System::compute_energy() const {
    std::vector<std::vector<std::array<double, 3>>> forces;
    forces.resize(residues_.size());
    for (size_t i = 0; i < residues_.size(); ++i) {
        forces[i].resize(residues_[i].particles.size(), {0.0, 0.0, 0.0});
    }
    
    double total_energy = 0.0;
    for (const auto& force : forces_) {
        total_energy += force->calculate_forces(*this, forces);
    }
    return total_energy;
}

void System::update_positions(double dt) {
    for (auto& residue : residues_) {
        for (auto& particle : residue.particles) {
            if (!particle.is_virtual) {
                for (size_t i = 0; i < 3; ++i) {
                    particle.position[i] += particle.velocity[i] * dt;
                }
            }
        }
    }
}

void System::update_velocities(double dt) {
    std::vector<std::vector<std::array<double, 3>>> forces;
    forces.resize(residues_.size());
    for (size_t i = 0; i < residues_.size(); ++i) {
        forces[i].resize(residues_[i].particles.size(), {0.0, 0.0, 0.0});
    }
    
    // Calculate forces
    for (const auto& force : forces_) {
        force->calculate_forces(*this, forces);
    }
    
    // Update velocities
    for (size_t i = 0; i < residues_.size(); ++i) {
        for (size_t j = 0; j < residues_[i].particles.size(); ++j) {
            auto& particle = residues_[i].particles[j];
            if (!particle.is_virtual && particle.mass > 0.0) {
                for (size_t k = 0; k < 3; ++k) {
                    particle.velocity[k] += forces[i][j][k] * dt / particle.mass;
                }
            }
        }
    }
}

// Distance computation
std::array<double, 3> System::compute_distance(const Particle& p1, const Particle& p2) const {
    std::array<double, 3> dr;
    for (size_t i = 0; i < 3; ++i) {
        dr[i] = p2.position[i] - p1.position[i];
    }

    if (has_periodic_boundary_) {
        // Convert to fractional coordinates
        std::array<double, 3> s = {0.0, 0.0, 0.0};
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                s[i] += dr[j] * box_vectors_[i][j];
            }
        }

        // Apply minimum image convention
        for (double& x : s) {
            x -= std::round(x);
        }

        // Convert back to Cartesian coordinates
        std::fill(dr.begin(), dr.end(), 0.0);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                dr[i] += s[j] * box_vectors_[j][i];
            }
        }
    }

    return dr;
}

// Private validation methods
void System::validate_residue_index(size_t index) const {
    if (index >= residues_.size()) {
        throw std::out_of_range("Residue index out of range");
    }
}

void System::validate_particle_index(size_t residue_index, size_t particle_index) const {
    validate_residue_index(residue_index);
    if (particle_index >= residues_[residue_index].particles.size()) {
        throw std::out_of_range("Particle index out of range");
    }
}

void System::validate_constraint_index(size_t index) const {
    if (index >= constraints_.size()) {
        throw std::out_of_range("Constraint index out of range");
    }
}

void System::validate_force_index(size_t index) const {
    if (index >= forces_.size()) {
        throw std::out_of_range("Force index out of range");
    }
}

void System::validate_box_vectors(const std::array<double, 3>& a,
                                const std::array<double, 3>& b,
                                const std::array<double, 3>& c) const {
    // Check that box vectors are finite
    for (const auto& v : {a, b, c}) {
        for (double x : v) {
            if (!std::isfinite(x)) {
                throw std::invalid_argument("Box vectors must be finite");
            }
        }
    }
    
    // Check that box vectors form a valid triclinic box
    double volume = a[0] * (b[1] * c[2] - b[2] * c[1]) -
                   a[1] * (b[0] * c[2] - b[2] * c[0]) +
                   a[2] * (b[0] * c[1] - b[1] * c[0]);
    if (volume <= 0.0) {
        throw std::invalid_argument("Box vectors must form a valid triclinic box with positive volume");
    }
}

void System::load_structure(const Structure& structure) {
    // Clear existing data
    residues_.clear();
    constraints_.clear();
    forces_.clear();
    has_periodic_boundary_ = false;

    // Get box information
    auto box = structure.get_box();
    if (box.has_value()) {
        auto box_vectors = structure.get_box_vectors();
        set_periodic_box_vectors(box_vectors[0], box_vectors[1], box_vectors[2]);
    }
    
    // Get atoms data and transfer to system
    auto atoms_data = structure.get_atoms_data();
    
    // Group atoms by residue
    std::map<std::pair<std::string, int>, std::vector<std::unordered_map<std::string, std::variant<std::string, int, double>>>> residue_atoms;
    for (const auto& atom : atoms_data) {
        std::string residue_name = std::get<std::string>(atom.at("residue"));
        int sequence = std::get<int>(atom.at("sequence"));
        residue_atoms[{residue_name, sequence}].push_back(atom);
    }
    
    // Create residues and add atoms
    for (const auto& [residue_key, atoms] : residue_atoms) {
        const auto& [residue_name, sequence] = residue_key;
        size_t res_idx = add_residue(residue_name);
        
        for (const auto& atom : atoms) {
            Particle p;
            p.position = {
                std::get<double>(atom.at("x")),
                std::get<double>(atom.at("y")),
                std::get<double>(atom.at("z"))
            };
            p.velocity = {0.0, 0.0, 0.0};  // Initialize velocities to zero
            
            // Try to get mass and charge from different possible fields
            bool mass_found = false;
            std::vector<std::string> mass_fields = {"mass", "topo_mass", "atom_mass"};
            for (const auto& field : mass_fields) {
                try {
                    if (atom.find(field) != atom.end()) {
                        p.mass = std::get<double>(atom.at(field));
                        if (std::isfinite(p.mass) && p.mass > 0.0) {
                            mass_found = true;
                            break;
                        }
                    }
                } catch (const std::exception&) {
                    continue;
                }
            }

            // Try to get charge from different possible fields
            std::vector<std::string> charge_fields = {"charge", "topo_charge", "atom_charge"};
            for (const auto& field : charge_fields) {
                try {
                    if (atom.find(field) != atom.end()) {
                        p.charge = std::get<double>(atom.at(field));
                        if (std::isfinite(p.charge)) {
                            break;
                        }
                    }
                } catch (const std::exception&) {
                    continue;
                }
            }

            // If no valid mass found, use a default mass
            if (!mass_found) {
                p.mass = 1.0;  // Default mass in atomic mass units
            }
            
            p.is_virtual = false;  // Default to non-virtual
            add_particle(res_idx, p);
        }
    }
}

void System::load_structure_psf(const std::string& pdb_file, const std::string& psf_file) {
    // Create a temporary project to load the structure
    Project project("temp_project");
    auto structure = project.create_structure();

    // Load PDB file
    structure.read_pdb(pdb_file);

    // Try standard PSF loading first
    try {
        structure.read_psf(psf_file);
    } catch (const std::runtime_error& e) {
        // If standard method fails, try using multiple PSF method
        std::vector<io::PDBAtom*> atom_ptrs;
        for (const auto& residue : structure.residues()) {
            for (const auto& atom : residue->atom_ptrs) {
                if (atom) {
                    atom_ptrs.push_back(atom.get());
                }
            }
        }

        // Try multi-residue approach first
        int updated = io::PSFParser::update_pdb_atoms_multi_residue(atom_ptrs, psf_file);
        if (updated == 0) {
            // If multi-residue approach fails, try single-residue approach
            updated = io::PSFParser::update_pdb_atoms_single_residue(atom_ptrs, psf_file, "BENX");
        }
        
        if (updated == 0) {
            throw std::runtime_error("No atoms were updated with PSF information using any method");
        }
    }

    // Load the structure into the system
    load_structure(structure);
}

void System::load_structure_psf_multi(const std::string& pdb_file, const std::string& psf_file) {
    // Create a temporary project to load the structure
    Project project("temp_project");
    auto structure = project.create_structure();

    // Load PDB file
    structure.read_pdb(pdb_file);

    // Get all atom pointers
    std::vector<io::PDBAtom*> atom_ptrs;
    for (const auto& residue : structure.residues()) {
        for (const auto& atom : residue->atom_ptrs) {
            if (atom) {
                atom_ptrs.push_back(atom.get());
            }
        }
    }

    // Try to update atoms using multi-residue method
    int updated = io::PSFParser::update_pdb_atoms_multi_residue(atom_ptrs, psf_file);
    if (updated == 0) {
        throw std::runtime_error("No atoms were updated with PSF information using multi-residue method");
    }

    // Load the structure into the system
    load_structure(structure);
}

void System::load_structure_psf_single(const std::string& pdb_file, const std::string& psf_file, const std::string& target_residue) {
    // Create a temporary project to load the structure
    Project project("temp_project");
    auto structure = project.create_structure();

    // Load PDB file
    structure.read_pdb(pdb_file);

    // Get all atom pointers
    std::vector<io::PDBAtom*> atom_ptrs;
    for (const auto& residue : structure.residues()) {
        for (const auto& atom : residue->atom_ptrs) {
            if (atom) {
                atom_ptrs.push_back(atom.get());
            }
        }
    }

    // Try to update atoms using single-residue method
    int updated = io::PSFParser::update_pdb_atoms_single_residue(atom_ptrs, psf_file, target_residue);
    if (updated == 0) {
        throw std::runtime_error("No atoms were updated with PSF information using single-residue method");
    }

    // Load the structure into the system
    load_structure(structure);
}

void System::load_structure_top(const std::string& pdb_file, const std::string& top_file) {
    // Create a temporary project to load the structure
    Project project("temp_project");
    auto structure = project.create_structure();

    // Load PDB and TOP files
    structure.read_pdb(pdb_file);
    structure.read_top(top_file);

    // Load the structure into the system
    load_structure(structure);
}

void System::load_structure_psf_auto(const std::string& pdb_file, const std::string& psf_file) {
    // Create a temporary project to load the structure
    Project project("temp_project");
    auto structure = project.create_structure();

    // Load PDB file
    structure.read_pdb(pdb_file);

    // Get all atom pointers
    std::vector<io::PDBAtom*> atom_ptrs;
    for (const auto& residue : structure.residues()) {
        for (const auto& atom : residue->atom_ptrs) {
            if (atom) {
                atom_ptrs.push_back(atom.get());
            }
        }
    }

    // Try multi-residue approach first
    int updated = io::PSFParser::update_pdb_atoms_multi_residue(atom_ptrs, psf_file);
    if (updated > 0) {
        // Multi-residue approach worked, load the structure
        load_structure(structure);
        return;
    }

    // If multi-residue approach failed, try single-residue approach with different residue names
    std::set<std::string> residue_names;
    for (const auto& residue : structure.residues()) {
        residue_names.insert(residue->name);
    }

    for (const auto& residue_name : residue_names) {
        updated = io::PSFParser::update_pdb_atoms_single_residue(atom_ptrs, psf_file, residue_name);
        if (updated > 0) {
            // Single-residue approach worked with this residue name
            load_structure(structure);
            return;
        }
    }

    // If both approaches failed, throw an error
    throw std::runtime_error("Failed to load PSF file: neither multi-residue nor single-residue approach worked");
}

} // namespace core
} // namespace pygcmc
