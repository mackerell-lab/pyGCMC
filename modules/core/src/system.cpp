// modules/core/src/system.cpp

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

namespace {
// Single, unified helper function to convert IO::PDBAtom to core::Particle
Particle convert_io_to_core_particle(const io::PDBAtom& io_atom) {
    Particle particle;
    particle.serial = io_atom.serial;
    particle.name = io_atom.name;
    particle.residue = io_atom.residue;
    particle.sequence = io_atom.sequence;
    particle.x = io_atom.x;
    particle.y = io_atom.y;
    particle.z = io_atom.z;
    particle.charge = std::stod(io_atom.charge); // Ensure charge is converted correctly
    particle.type = io_atom.type;
    particle.nameTop = io_atom.element;
    return particle;
}

// Helper function to convert IO::IOResidue to core::Residue
Residue convert_io_to_core_residue(const io::IOResidue& io_residue) {
    Residue core_residue;
    core_residue.name = io_residue.name;
    core_residue.sequence_number = io_residue.sequence_number;
    core_residue.chain_id = io_residue.chain_id;
    
    core_residue.atoms.reserve(io_residue.atoms.size());
    for (const auto& io_atom : io_residue.atoms) {
        core_residue.atoms.push_back(convert_io_to_core_particle(io_atom));
    }
    
    return core_residue;
}
}  // anonymous namespace

// File loading methods
void System::load_pdb(const std::string& filename) {
    try {
        auto parsed = io::PDBParser::parse(filename);
        auto cryst = parsed.first;
        auto parsed_residues = parsed.second;

        std::cout << "Successfully loaded " << parsed_residues.size() 
                  << " residues from " << filename << std::endl;
        std::cout << "Cell parameters: a=" << cryst[0] << ", b=" 
                  << cryst[1] << ", c=" << cryst[2] << std::endl;

        for (const auto& io_residue : parsed_residues) {
            add_residue(convert_io_to_core_residue(io_residue));
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
        for (size_t i = 0; i < residues_.size(); ++i) {
            for (size_t j = 0; j < residues_[i].atoms.size() && j < top.atom_types.size(); ++j) {
                residues_[i].atoms[j].type = top.atom_types[j].name;
                residues_[i].atoms[j].charge = top.atom_types[j].charge;
                residues_[i].atoms[j].nameTop = top.atom_types[j].name;
            }
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
        for (size_t i = 0; i < residues_.size(); ++i) {
            for (size_t j = 0; j < residues_[i].atoms.size() && j < itp_atoms.size(); ++j) {
                residues_[i].atoms[j].type = itp_atoms[j].type;
                residues_[i].atoms[j].charge = itp_atoms[j].charge;
                residues_[i].atoms[j].nameTop = itp_atoms[j].name;
            }
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

// Residue management
void System::add_residue(const Residue& residue) {
    if (!residue.is_valid()) {
        throw SystemError("Invalid residue cannot be added to the system.");
    }
    residues_.emplace_back(residue);
}

void System::remove_residue(int index) {
    if (index >= 0 && index < static_cast<int>(residues_.size())) {
        residues_.erase(residues_.begin() + index);
    } else {
        throw std::out_of_range("Residue index out of range");
    }
}

size_t System::get_residue_count() const {
    return residues_.size();
}

const Residue& System::get_residue(size_t index) const {
    if (index >= residues_.size()) {
        throw std::out_of_range("Residue index out of range");
    }
    return residues_[index];
}

Residue& System::get_residue(size_t index) {
    if (index >= residues_.size()) {
        throw std::out_of_range("Residue index out of range");
    }
    return residues_[index];
}

// Energy computation
double System::compute_total_energy() const {
    double total_energy = 0.0;
    
    for (size_t i = 0; i < residues_.size(); ++i) {
        for (size_t j = i + 1; j < residues_.size(); ++j) {
            total_energy += compute_residue_energy(residues_[i], residues_[j]);
        }
    }
    
    return total_energy;
}

double System::compute_residue_energy(const Residue& res1, const Residue& res2) const {
    double energy = 0.0;
    for (const auto& atom1 : res1.atoms) {
        for (const auto& atom2 : res2.atoms) {
            energy += compute_pair_energy(atom1, atom2);
        }
    }
    return energy;
}

std::pair<double, double> System::get_system_state() const {
    double kinetic_energy = 0.0;
    double potential_energy = compute_total_energy();
    
    // Calculate kinetic energy
    for (const auto& res : residues_) {
        for (const auto& p : res.atoms) {
            double v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
            kinetic_energy += 0.5 * v2;  // Assuming mass = 1 for simplicity
        }
    }
    
    return {kinetic_energy, potential_energy};
}

// Dynamics methods
void System::update_positions(double dt) {
    for (auto& res : residues_) {
        for (auto& particle : res.atoms) {
            particle.x += particle.vx * dt;
            particle.y += particle.vy * dt;
            particle.z += particle.vz * dt;

            if (use_periodic_) {
                particle.apply_periodic_boundary(box_size_);
            }
        }
    }
}

void System::update_velocities(double dt) {
    std::vector<std::array<double, 3>> accelerations(residues_.size(), {0.0, 0.0, 0.0});
    
    for (size_t i = 0; i < residues_.size(); ++i) {
        for (size_t j = i + 1; j < residues_.size(); ++j) {
            for (const auto& atom1 : residues_[i].atoms) {
                for (const auto& atom2 : residues_[j].atoms) {
                    double dx = atom1.x - atom2.x;
                    double dy = atom1.y - atom2.y;
                    double dz = atom1.z - atom2.z;
                    double distance = std::sqrt(dx * dx + dy * dy + dz * dz) + 1e-12;

                    double inv_r = sigma_ / distance;
                    double inv_r6 = std::pow(inv_r, 6);
                    double inv_r12 = std::pow(inv_r6, 2);

                    // Calculate forces
                    double force_lj = 24.0 * epsilon_ / distance * (2.0 * inv_r12 - inv_r6);
                    double force_coulomb = (atom1.charge * atom2.charge) / (distance * distance);
                    double total_force = force_lj + force_coulomb;

                    // Update accelerations
                    double fx = total_force * (dx / distance);
                    double fy = total_force * (dy / distance);
                    double fz = total_force * (dz / distance);

                    // 累加力
                    // 由于加速度是基于残基的，这里需要按残基索引累加
                    accelerations[i][0] += fx;
                    accelerations[i][1] += fy;
                    accelerations[i][2] += fz;

                    accelerations[j][0] -= fx;
                    accelerations[j][1] -= fy;
                    accelerations[j][2] -= fz;
                }
            }
        }
    }

    // 更新速度
    for (size_t i = 0; i < residues_.size(); ++i) {
        for (auto& atom : residues_[i].atoms) {
            atom.vx += accelerations[i][0] * dt;
            atom.vy += accelerations[i][1] * dt;
            atom.vz += accelerations[i][2] * dt;
        }
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
    return x - box_size_ * std::floor(x / box_size_);
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
