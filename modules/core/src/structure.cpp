// modules/core/src/structure.cpp

#include "pygcmc/core/structure.hpp"
#include "pygcmc/core/forcefield.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/io/itp_parser.hpp"
#include <memory>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <map>
#include <algorithm>
#include <sstream>

namespace pygcmc {
namespace core {

void Structure::apply_forcefield(const std::shared_ptr<ForceField>& forcefield) {
    if (!forcefield) {
        throw std::runtime_error("Cannot apply null forcefield");
    }
    
    // Store the forcefield
    forcefield_ = forcefield;
    
    // Apply forcefield parameters to each atom
    for (auto& atom : atoms_) {
        if (atom && !atom->topo_type.empty()) {
            // Get parameters for this atom type
            const auto& params = forcefield->nonbonded_params();
            auto it = params.find(atom->topo_type);
            if (it != params.end()) {
                atom->forcefield_epsilon = it->second.epsilon;
                atom->forcefield_rmin = it->second.rmin;
            }
        }
    }
}

void Structure::add_residue(std::shared_ptr<io::IOResidue> residue) {
    if (residue) {
        residues_.push_back(residue);
    }
}

void Structure::add_atom(std::shared_ptr<io::PDBAtom> atom) {
    if (atom) {
        atoms_.push_back(atom);
    }
}

std::vector<std::array<double, 3>> Structure::get_coordinates() const {
    std::vector<std::array<double, 3>> coords;
    coords.reserve(atoms_.size());
    
    for (const auto& atom : atoms_) {
        if (atom) {
            coords.push_back({atom->x, atom->y, atom->z});
        }
    }
    
    return coords;
}

std::array<std::array<double, 3>, 3> Structure::get_box_vectors() const {
    std::array<std::array<double, 3>, 3> vectors{};
    
    if (box_) {
        const auto& [a, b, c, alpha, beta, gamma] = box_.value();
        
        // Convert angles from degrees to radians
        const double alpha_rad = alpha * M_PI / 180.0;
        const double beta_rad = beta * M_PI / 180.0;
        const double gamma_rad = gamma * M_PI / 180.0;
        
        // Calculate box vectors following the triclinic box convention
        vectors[0] = {a, 0.0, 0.0};
        vectors[1] = {b * cos(gamma_rad), b * sin(gamma_rad), 0.0};
        
        const double cx = c * cos(beta_rad);
        const double cy = c * (cos(alpha_rad) - cos(beta_rad) * cos(gamma_rad)) / sin(gamma_rad);
        const double cz = sqrt(c * c - cx * cx - cy * cy);
        vectors[2] = {cx, cy, cz};
    }
    
    return vectors;
}

std::unordered_map<std::string, double> Structure::get_energy_components() const {
    std::unordered_map<std::string, double> components;
    
    // Initialize all required components
    components["bond"] = 0.0;
    components["angle"] = 0.0;
    components["dihedral"] = 0.0;
    components["improper"] = 0.0;
    components["vdw"] = 0.0;
    components["electrostatic"] = 0.0;
    
    // Calculate energies if force field is available
    if (forcefield_) {
        // TODO: Implement actual energy calculations
        // For now, return placeholder values
        components["bond"] = 100.0;  // Example value
        components["angle"] = 200.0;  // Example value
        components["dihedral"] = 150.0;  // Example value
        components["improper"] = 50.0;  // Example value
        components["vdw"] = 300.0;  // Example value
        components["electrostatic"] = 400.0;  // Example value
    }
    
    // Calculate total energy
    double total = 0.0;
    for (const auto& [component, energy] : components) {
        if (component != "total") {
            total += energy;
        }
    }
    components["total"] = total;
    
    return components;
}

std::vector<std::tuple<size_t, double, double, double>> Structure::get_atom_energy_contributions() const {
    std::vector<std::tuple<size_t, double, double, double>> atom_energies;
    atom_energies.reserve(atoms_.size());
    
    if (!forcefield_) {
        return atom_energies;
    }
    
    // Calculate per-atom energy contributions
    for (size_t i = 0; i < atoms_.size(); ++i) {
        double vdw_energy = 0.0;
        double elec_energy = 0.0;
        
        const auto& atom1 = atoms_[i];
        if (!atom1) continue;
        
        for (size_t j = 0; j < atoms_.size(); ++j) {
            if (i == j) continue;
            
            const auto& atom2 = atoms_[j];
            if (!atom2) continue;
            
            // Calculate distance
            const double dx = atom2->x - atom1->x;
            const double dy = atom2->y - atom1->y;
            const double dz = atom2->z - atom1->z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r = sqrt(r2);
            
            // VDW energy contribution
            if (!std::isnan(atom1->forcefield_epsilon) && !std::isnan(atom2->forcefield_epsilon) &&
                !std::isnan(atom1->forcefield_rmin) && !std::isnan(atom2->forcefield_rmin)) {
                const double eps = sqrt(atom1->forcefield_epsilon * atom2->forcefield_epsilon);
                const double rmin = atom1->forcefield_rmin + atom2->forcefield_rmin;
                const double ratio = rmin / r;
                const double ratio6 = ratio * ratio * ratio * ratio * ratio * ratio;
                vdw_energy += 0.5 * eps * (ratio6 * ratio6 - 2.0 * ratio6);  // Half because of double counting
            }
            
            // Electrostatic energy contribution
            if (!std::isnan(atom1->topo_charge) && !std::isnan(atom2->topo_charge)) {
                const double k_coulomb = 332.0716;  // kcal/mol/Å/e^2
                elec_energy += 0.5 * k_coulomb * atom1->topo_charge * atom2->topo_charge / r;  // Half because of double counting
            }
        }
        
        atom_energies.emplace_back(i, vdw_energy, elec_energy, vdw_energy + elec_energy);
    }
    
    return atom_energies;
}

void Structure::load_pdb(const std::string& pdb_file) {
    read_pdb_file(pdb_file);
}

void Structure::load_top(const std::string& top_file) {
    read_top_file(top_file);
}

void Structure::load_top_with_includes(const std::string& top_file) {
    read_top_file_with_includes(top_file);
}

void Structure::update_atoms_topology(io::TopParser& top_parser) {
    // Get all atom pointers
    std::vector<io::PDBAtom*> atom_ptrs;
    for (const auto& residue : residues_) {
        for (const auto& atom : residue->atom_ptrs) {
            if (atom) {
                atom_ptrs.push_back(atom.get());
            }
        }
    }
    
    // Update atoms with topology information
    int updated = top_parser.update_pdb_atoms(atom_ptrs);
    if (updated == 0) {
        throw std::runtime_error("No atoms were updated with topology information");
    }
}

void Structure::read_pdb_file(const std::string& pdb_file) {
    // Parse PDB file
    auto [box_vec, residues] = io::PDBParser::parse(pdb_file);
    
    // Convert vector box to array box if present
    if (box_vec && box_vec->size() == 6) {
        std::array<double, 6> box_arr;
        std::copy(box_vec->begin(), box_vec->end(), box_arr.begin());
        set_box(std::make_optional(box_arr));
    } else {
        set_box(std::nullopt);
    }
    
    // Store existing topology information
    std::vector<std::tuple<std::string, int, std::string, std::string, double, double>> cached_topology;
    for (const auto& atom : atoms_) {
        if (atom) {
            cached_topology.emplace_back(atom->residue, atom->sequence, atom->name, 
                                       atom->topo_type, atom->topo_charge, atom->topo_mass);
        }
    }
    
    // Clear existing data
    residues_.clear();
    atoms_.clear();
    
    // Add residues to structure
    for (const auto& residue : residues) {
        // Create a new shared_ptr to a copy of the residue
        auto residue_ptr = std::make_shared<io::IOResidue>();
        *residue_ptr = residue;  // Use copy assignment

        // Create shared_ptr for each atom and update atom_ptrs
        residue_ptr->atom_ptrs.clear();  // Clear existing pointers
        for (const auto& atom : residue.atoms) {
            auto atom_ptr = std::make_shared<io::PDBAtom>(atom);
            
            // Restore topology information if it exists
            for (const auto& [res, seq, name, type, charge, mass] : cached_topology) {
                if (res == atom.residue && seq == atom.sequence && name == atom.name) {
                    atom_ptr->topo_type = type;
                    atom_ptr->topo_charge = charge;
                    atom_ptr->topo_mass = mass;
                    break;
                }
            }
            
            residue_ptr->atom_ptrs.push_back(atom_ptr);
            add_atom(atom_ptr);  // Add atom to structure's atoms_ vector
        }
        
        add_residue(residue_ptr);
    }

    // Apply cached topology in order: other topology first, then ITPs, and PSF last to take precedence
    apply_cached_topology();
    apply_cached_itps();
    if (has_cached_psf_ && cached_psf_) {
        apply_cached_psf();
    }
}

void Structure::read_top_file(const std::string& top_file) {
    // Default behavior is to read with includes
    read_top_file_with_includes(top_file);
}

void Structure::read_top_file_without_includes(const std::string& top_file) {
    // Create TopParser instance and parse the file without includes
    io::TopParser top_parser;
    if (!top_parser.parse(top_file)) {
        throw std::runtime_error("Failed to parse topology file: " + top_file);
    }
    
    if (atoms_.empty()) {
        // Cache the topology for later use
        cached_topology_ = std::move(top_parser);
        has_cached_topology_ = true;
    } else {
        // Update atoms with topology information
        update_atoms_topology(top_parser);
    }
}

void Structure::read_top_file_with_includes(const std::string& top_file) {
    // Create TopParser instance and parse the file with includes
    io::TopParser top_parser;
    if (!top_parser.parse_with_includes(top_file)) {
        throw std::runtime_error("Failed to parse topology file: " + top_file);
    }
    
    if (atoms_.empty()) {
        // Cache the topology for later use
        cached_topology_ = std::move(top_parser);
        has_cached_topology_ = true;
    } else {
        // Update atoms with topology information
        update_atoms_topology(top_parser);
    }
}

void Structure::apply_cached_topology() {
    if (has_cached_topology_ && !atoms_.empty()) {
        update_atoms_topology(*cached_topology_);
        cached_topology_ = std::nullopt;
        has_cached_topology_ = false;
    }
}

std::vector<std::unordered_map<std::string, std::variant<std::string, int, double>>> Structure::get_atoms_data() const {
    std::vector<std::unordered_map<std::string, std::variant<std::string, int, double>>> atoms_data;
    atoms_data.reserve(atoms_.size());

    for (const auto& atom : atoms_) {
        if (!atom) continue;

        std::unordered_map<std::string, std::variant<std::string, int, double>> atom_data;
        
        // Basic atom information
        atom_data["name"] = atom->name;
        atom_data["residue"] = atom->residue;
        atom_data["sequence"] = atom->sequence;
        atom_data["type"] = atom->type;
        
        // Coordinates
        atom_data["x"] = atom->x;
        atom_data["y"] = atom->y;
        atom_data["z"] = atom->z;
        
        // Topology information
        atom_data["topo_type"] = atom->topo_type;
        atom_data["topo_charge"] = atom->topo_charge;
        atom_data["topo_mass"] = atom->topo_mass;
        
        // Force field parameters if available
        if (!std::isnan(atom->forcefield_epsilon)) {
            atom_data["forcefield_epsilon"] = atom->forcefield_epsilon;
        }
        if (!std::isnan(atom->forcefield_rmin)) {
            atom_data["forcefield_rmin"] = atom->forcefield_rmin;
        }
        
        atoms_data.push_back(std::move(atom_data));
    }
    
    return atoms_data;
}

void Structure::load_psf(const std::string& psf_file) {
    read_psf_file(psf_file);
}

void Structure::load_itp(const std::string& itp_file) {
    read_itp_file(itp_file);
}

void Structure::read_psf_file(const std::string& psf_file) {
    // Create PSFParser instance and parse the file
    io::PSFParser psf_parser;
    try {
        if (!psf_parser.parse(psf_file)) {
            throw std::runtime_error("Failed to parse PSF file: " + psf_file);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse PSF file: " + psf_file + " - " + e.what());
    }
    
    if (atoms_.empty()) {
        // Cache the PSF for later use
        cached_psf_ = std::move(psf_parser);
        has_cached_psf_ = true;
    } else {
        // Get all atom pointers
        std::vector<io::PDBAtom*> atom_ptrs;
        for (const auto& residue : residues_) {
            for (const auto& atom : residue->atom_ptrs) {
                if (atom) {
                    atom_ptrs.push_back(atom.get());
                }
            }
        }
        
        // First try standard method (order-based mapping)
        int updated = psf_parser.update_pdb_atoms_by_order(atom_ptrs);
        if (updated == 0) {
            // If standard method fails, try updating with multiple PSF method
            std::vector<std::string> psf_files = {psf_file};
            updated = io::PSFParser::update_pdb_atoms_from_multiple_psf(atom_ptrs, psf_files);
            
            if (updated == 0) {
                throw std::runtime_error("No atoms were updated with PSF information using either method");
            }
        }
    }
}

void Structure::read_itp_file(const std::string& itp_file) {
    // Create ITPParser instance and parse the file
    io::ITPParser itp_parser;
    if (!itp_parser.parse(itp_file)) {
        throw std::runtime_error("Cannot open file: " + itp_file);
    }
    
    if (atoms_.empty()) {
        // Cache the ITP for later use
        cached_itps_.push_back(std::move(itp_parser));
    } else {
        // Update atoms with ITP information
        update_atoms_topology(itp_parser);
    }
}

void Structure::update_atoms_topology(io::PSFParser& psf_parser) {
    // Get all atom pointers
    std::vector<io::PDBAtom*> atom_ptrs;
    for (const auto& residue : residues_) {
        for (const auto& atom : residue->atom_ptrs) {
            if (atom) {
                atom_ptrs.push_back(atom.get());
            }
        }
    }
    
    // Update atoms with PSF information
    int updated = psf_parser.update_pdb_atoms(atom_ptrs);
    if (updated == 0) {
        throw std::runtime_error("No atoms were updated with PSF information");
    }
}

void Structure::update_atoms_topology(io::ITPParser& itp_parser) {
    // Get all atom pointers
    std::vector<io::PDBAtom*> atom_ptrs;
    for (const auto& residue : residues_) {
        for (const auto& atom : residue->atom_ptrs) {
            if (atom) {
                atom_ptrs.push_back(atom.get());
            }
        }
    }
    
    // Update atoms with ITP information
    // Don't check number of updated atoms - this is expected for ITP files
    // that only contain some residues
    itp_parser.update_pdb_atoms(atom_ptrs);
}

void Structure::apply_cached_psf() {
    if (has_cached_psf_ && !atoms_.empty() && cached_psf_) {
        update_atoms_topology(*cached_psf_);
        cached_psf_ = std::nullopt;
        has_cached_psf_ = false;
    }
}

void Structure::apply_cached_itps() {
    if (!atoms_.empty()) {
        for (auto& itp : cached_itps_) {
            update_atoms_topology(itp);
        }
        cached_itps_.clear();
    }
}

void Structure::load_structure_psf_single(const std::string& psf_file, const std::string& pdb_file) {
    // Load PDB first to get coordinates
    load_structure_pdb(pdb_file);
    
    // Now parse PSF file to get topology information
    std::ifstream psf(psf_file);
    if (!psf.is_open()) {
        throw std::runtime_error("Could not open PSF file: " + psf_file);
    }
    
    std::string line;
    bool reading_atoms = false;
    int num_atoms = 0;
    int atom_count = 0;
    
    while (std::getline(psf, line)) {
        if (line.find("!NATOM") != std::string::npos) {
            std::istringstream iss(line);
            iss >> num_atoms;
            reading_atoms = true;
            continue;
        }
        
        if (reading_atoms && atom_count < num_atoms) {
            // PSF format: ID SEGID RESID RES ATOM TYPE CHARGE MASS 0
            std::istringstream iss(line);
            int id;
            std::string segid, resid, res, atom, type;
            double charge, mass;
            int zero;
            
            iss >> id >> segid >> resid >> res >> atom >> type >> charge >> mass >> zero;
            
            // Find corresponding atom in atoms_data_
            for (auto& atom_data : atoms_data_) {
                if (std::get<std::string>(atom_data["residue"]) == res &&
                    std::get<int>(atom_data["sequence"]) == std::stoi(resid) &&
                    std::get<std::string>(atom_data["name"]) == atom) {
                    // Store topology information
                    atom_data["topo_mass"] = mass;
                    atom_data["topo_charge"] = charge;
                    atom_data["topo_type"] = type;
                    break;
                }
            }
            
            atom_count++;
        }
        
        if (atom_count >= num_atoms) {
            reading_atoms = false;
        }
    }
}

void Structure::load_structure_psf_auto(const std::string& psf_file, const std::string& pdb_file) {
    // Load PDB first to get coordinates
    load_structure_pdb(pdb_file);
    
    // Now parse PSF file to get topology information
    std::ifstream psf(psf_file);
    if (!psf.is_open()) {
        throw std::runtime_error("Could not open PSF file: " + psf_file);
    }
    
    std::string line;
    bool reading_atoms = false;
    int num_atoms = 0;
    int atom_count = 0;
    
    // Map to store topology information by residue name and atom name
    std::map<std::pair<std::string, std::string>, std::tuple<double, double, std::string>> topo_info;
    
    while (std::getline(psf, line)) {
        if (line.find("!NATOM") != std::string::npos) {
            std::istringstream iss(line);
            iss >> num_atoms;
            reading_atoms = true;
            continue;
        }
        
        if (reading_atoms && atom_count < num_atoms) {
            // PSF format: ID SEGID RESID RES ATOM TYPE CHARGE MASS 0
            std::istringstream iss(line);
            int id;
            std::string segid, resid, res, atom, type;
            double charge, mass;
            int zero;
            
            iss >> id >> segid >> resid >> res >> atom >> type >> charge >> mass >> zero;
            
            // Store topology information by residue and atom name
            topo_info[{res, atom}] = std::make_tuple(mass, charge, type);
            
            atom_count++;
        }
        
        if (atom_count >= num_atoms) {
            reading_atoms = false;
        }
    }
    
    // Apply topology information to all matching atoms
    for (auto& atom_data : atoms_data_) {
        std::string res = std::get<std::string>(atom_data["residue"]);
        std::string atom = std::get<std::string>(atom_data["name"]);
        
        auto it = topo_info.find({res, atom});
        if (it != topo_info.end()) {
            const auto& [mass, charge, type] = it->second;
            atom_data["topo_mass"] = mass;
            atom_data["topo_charge"] = charge;
            atom_data["topo_type"] = type;
        }
    }
}

} // namespace core
} // namespace pygcmc 