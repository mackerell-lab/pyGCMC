// modules/core/src/system.cpp

#include "pygcmc/core/system.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>

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

} // namespace core
} // namespace pygcmc
