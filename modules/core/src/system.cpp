#include "pygcmc/core/system.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/itp_parser.hpp"
#include "pygcmc/core/io/ff_parser.hpp"
#include <iostream>
#include <array>
#include <cmath>

namespace pygcmc {
namespace core {

System::System(double epsilon, double sigma)
    : epsilon_(epsilon), sigma_(sigma) {
    if (epsilon < 0.0 || sigma < 0.0) {
        throw std::invalid_argument("Force field parameters must be non-negative");
    }
}

System::~System() = default;

// File loading methods
void System::load_pdb(const std::string& filename) {
    try {
        auto parsed = io::PDBParser::parse(filename);
        auto cryst = parsed.first;
        auto parsed_atoms = parsed.second;

        std::cout << "Successfully loaded " << parsed_atoms.size() 
                  << " particles from " << filename << std::endl;
        std::cout << "Cell parameters: a=" << cryst[0] << ", b=" 
                  << cryst[1] << ", c=" << cryst[2] << std::endl;

        for (const auto& p : parsed_atoms) {
            Particle particle;
            particle.serial = p.serial;
            particle.name = p.name;
            particle.residue = p.residue;
            particle.sequence = p.sequence;
            particle.x = p.x;
            particle.y = p.y;
            particle.z = p.z;
            particle.charge = p.charge;
            particle.type = p.type;
            particle.nameTop = "";
            particle.typeNum = 0;
            particle.vx = 0.0;
            particle.vy = 0.0;
            particle.vz = 0.0;
            add_particle(particle);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to load PDB file: " << e.what() << std::endl;
        throw;
    }
}

void System::load_psf(const std::string& filename) {
    try {
        auto psf = io::PSFParser::parse(filename);
        // 处理PSF中的键、角度等信息（当前仅处理键）
        std::cout << "Successfully loaded PSF file: " << filename << std::endl;
        std::cout << "Number of bonds: " << psf.bonds.size() << std::endl;

        // 可以将键信息存储在系统中，供能量计算使用
        // 这里省略具体实现
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to load PSF file: " << e.what() << std::endl;
    }
}

void System::load_top(const std::string& filename) {
    try {
        auto top = io::TopParser::parse(filename);
        std::cout << "Successfully loaded TOP file: " << filename << std::endl;
        std::cout << "Number of atom types: " << top.atom_types.size() << std::endl;

        // 将TOP中的原子类型信息与粒子关联
        for (size_t i = 0; i < particles_.size() && i < top.atom_types.size(); ++i) {
            particles_[i].type = top.atom_types[i].name;
            particles_[i].charge = top.atom_types[i].charge;
            particles_[i].nameTop = top.atom_types[i].name;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to load TOP file: " << e.what() << std::endl;
    }
}

void System::load_itp(const std::string& filename) {
    try {
        auto itp_atoms = io::ITPParser::parse(filename);
        std::cout << "Successfully loaded ITP file: " << filename << std::endl;
        std::cout << "Number of ITP atoms: " << itp_atoms.size() << std::endl;

        // 将ITP中的原子类型信息与粒子关联
        for (size_t i = 0; i < particles_.size() && i < itp_atoms.size(); ++i) {
            particles_[i].type = itp_atoms[i].type;
            particles_[i].charge = itp_atoms[i].charge;
            particles_[i].nameTop = itp_atoms[i].name;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to load ITP file: " << e.what() << std::endl;
    }
}

void System::load_forcefield(const std::string& filename) {
    try {
        auto ff = io::FFParser::parse(filename);
        nb_dict_ = ff.first;
        nbfix_dict_ = ff.second;
        std::cout << "Successfully loaded force field file: " << filename << std::endl;
        std::cout << "Number of non-bonded parameters: " << nb_dict_.size() 
                  << ", number of NBFIX parameters: " << nbfix_dict_.size() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to load force field file: " << e.what() << std::endl;
    }
}

// Particle management
void System::add_particle(const Particle& particle) {
    particles_.emplace_back(particle);
}

void System::remove_particle(int index) {
    if (index >= 0 && index < static_cast<int>(particles_.size())) {
        particles_.erase(particles_.begin() + index);
    }
}

size_t System::get_particle_count() const {
    return particles_.size();
}

const Particle& System::get_particle(size_t index) const {
    if (index >= particles_.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    return particles_[index];
}

Particle& System::get_particle(size_t index) {
    if (index >= particles_.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    return particles_[index];
}

// Energy computation
double System::compute_total_energy() const {
    double total_energy = 0.0;
    
    for (size_t i = 0; i < particles_.size(); ++i) {
        for (size_t j = i + 1; j < particles_.size(); ++j) {
            total_energy += compute_pair_energy(particles_[i], particles_[j]);
        }
    }
    
    return total_energy;
}

std::pair<double, double> System::get_system_state() const {
    double kinetic_energy = 0.0;
    double potential_energy = compute_total_energy();
    
    // Calculate kinetic energy
    for (const auto& p : particles_) {
        double v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
        kinetic_energy += 0.5 * v2;  // Assuming mass = 1 for simplicity
    }
    
    return {kinetic_energy, potential_energy};
}

// Dynamics methods
void System::update_positions(double dt) {
    for (auto& particle : particles_) {
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
        particle.z += particle.vz * dt;

        if (use_periodic_) {
            particle.apply_periodic_boundary(box_size_);
        }
    }
}

void System::update_velocities(double dt) {
    std::vector<std::array<double, 3>> accelerations(particles_.size(), {0.0, 0.0, 0.0});

    for (size_t i = 0; i < particles_.size(); ++i) {
        for (size_t j = i + 1; j < particles_.size(); ++j) {
            double dx = particles_[i].x - particles_[j].x;
            double dy = particles_[i].y - particles_[j].y;
            double dz = particles_[i].z - particles_[j].z;
            double distance = std::sqrt(dx * dx + dy * dy + dz * dz) + 1e-12;

            double inv_r = sigma_ / distance;
            double inv_r6 = std::pow(inv_r, 6);
            double inv_r12 = std::pow(inv_r6, 2);

            // Calculate forces
            double force_lj = 24.0 * epsilon_ / distance * (2.0 * inv_r12 - inv_r6);
            double force_coulomb = (particles_[i].charge * particles_[j].charge) / (distance * distance);
            double total_force = force_lj + force_coulomb;

            // Update accelerations
            double fx = total_force * (dx / distance);
            double fy = total_force * (dy / distance);
            double fz = total_force * (dz / distance);

            accelerations[i][0] += fx;
            accelerations[i][1] += fy;
            accelerations[i][2] += fz;

            accelerations[j][0] -= fx;
            accelerations[j][1] -= fy;
            accelerations[j][2] -= fz;
        }
    }

    // Update velocities
    for (size_t i = 0; i < particles_.size(); ++i) {
        particles_[i].vx += accelerations[i][0] * dt;
        particles_[i].vy += accelerations[i][1] * dt;
        particles_[i].vz += accelerations[i][2] * dt;
    }
}

// Boundary conditions
void System::set_periodic_boundary(double box_size) {
    if (box_size <= 0.0) {
        throw std::invalid_argument("Box size must be positive");
    }
    box_size_ = box_size;
    use_periodic_ = true;
}

double System::apply_pbc(double x) const {
    if (!use_periodic_) return x;
    return x - box_size_ * std::round(x / box_size_);
}

std::array<double, 3> System::compute_distance(const Particle& p1, const Particle& p2) const {
    std::array<double, 3> dr = {p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
    
    if (use_periodic_) {
        for (auto& d : dr) {
            d = apply_pbc(d);
        }
    }
    
    return dr;
}

double System::compute_pair_energy(const Particle& p1, const Particle& p2) const {
    auto dr = compute_distance(p1, p2);
    double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
    double distance = std::sqrt(r2);
    
    if (distance < 1e-10) return 0.0;  // Skip self-interaction or overlapping particles

    // Try to find NBFIX parameters first
    auto pair_key = std::make_pair(p1.nameTop, p2.nameTop);
    double lj_energy = 0.0;

    auto nbfix_it = nbfix_dict_.find(pair_key);
    if (nbfix_it != nbfix_dict_.end()) {
        // Use NBFIX parameters
        const auto& params = nbfix_it->second;
        double sigma = params.param1;
        double epsilon = params.param2;
        double inv_r = sigma / distance;
        double inv_r6 = std::pow(inv_r, 6);
        lj_energy = 4.0 * epsilon * (std::pow(inv_r6, 2) - inv_r6);
    } else {
        // Fall back to regular NB parameters with combining rules
        auto type1_it = nb_dict_.find(p1.nameTop);
        auto type2_it = nb_dict_.find(p2.nameTop);
        
        if (type1_it != nb_dict_.end() && type2_it != nb_dict_.end()) {
            // Use Lorentz-Berthelot combining rules
            double sigma = 0.5 * (type1_it->second.param1 + type2_it->second.param1);
            double epsilon = std::sqrt(type1_it->second.param2 * type2_it->second.param2);
            double inv_r = sigma / distance;
            double inv_r6 = std::pow(inv_r, 6);
            lj_energy = 4.0 * epsilon * (std::pow(inv_r6, 2) - inv_r6);
        } else {
            // Use default parameters
            double inv_r = sigma_ / distance;
            double inv_r6 = std::pow(inv_r, 6);
            lj_energy = 4.0 * epsilon_ * (std::pow(inv_r6, 2) - inv_r6);
        }
    }
    
    double coulomb_energy = (p1.charge * p2.charge) / distance;
    return lj_energy + coulomb_energy;
}

} // namespace core
} // namespace pygcmc
