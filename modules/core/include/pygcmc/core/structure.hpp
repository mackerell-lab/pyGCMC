// modules/core/include/pygcmc/core/structure.hpp

#pragma once

#include <memory>
#include <vector>
#include <optional>
#include <array>
#include <unordered_map>
#include <tuple>
#include <variant>
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/io/itp_parser.hpp"
#include "pygcmc/core/forcefield.hpp"

namespace pygcmc {
namespace core {

class ForceField;

class Structure {
public:
    Structure() = default;
    ~Structure() = default;

    // Box dimensions methods
    void set_box(const std::optional<std::array<double, 6>>& box) { box_ = box; }
    std::optional<std::array<double, 6>> get_box() const { return box_; }

    // Apply force field parameters to all atoms
    void apply_forcefield(const std::shared_ptr<ForceField>& forcefield);

    // Add a residue to the structure
    void add_residue(std::shared_ptr<io::IOResidue> residue);

    // Add an atom to the structure
    void add_atom(std::shared_ptr<io::PDBAtom> atom);

    // Get number of atoms
    size_t get_num_atoms() const { return atoms_.size(); }

    // Get residues
    const std::vector<std::shared_ptr<io::IOResidue>>& residues() const { return residues_; }

    // Get atoms
    const std::vector<std::shared_ptr<io::PDBAtom>>& atoms() const { return atoms_; }

    // New data access methods
    std::vector<std::array<double, 3>> get_coordinates() const;
    std::array<std::array<double, 3>, 3> get_box_vectors() const;
    
    // Energy calculation methods
    std::unordered_map<std::string, double> get_energy_components() const;
    std::vector<std::tuple<size_t, double, double, double>> get_atom_energy_contributions() const;

    // Get atoms data in a format suitable for Python
    std::vector<std::unordered_map<std::string, std::variant<std::string, int, double>>> get_atoms_data() const;

    // New functions for loading structure data
    void load_pdb(const std::string& pdb_file);
    void load_top(const std::string& top_file);
    void load_top_with_includes(const std::string& top_file);

    // New separate loading functions
    void read_pdb_file(const std::string& pdb_file);
    void read_top_file(const std::string& top_file);  // Default with includes
    void read_top_file_without_includes(const std::string& top_file);
    void read_top_file_with_includes(const std::string& top_file);  // Alias for read_top_file

    // Alias functions for test compatibility
    void read_pdb(const std::string& pdb_file) { read_pdb_file(pdb_file); }
    void read_top(const std::string& top_file) { read_top_file(top_file); }
    void read_top_without_includes(const std::string& top_file) { read_top_file_without_includes(top_file); }

    // New PSF and ITP loading functions
    void load_psf(const std::string& psf_file);
    void load_itp(const std::string& itp_file);
    
    // Alias methods for PSF and ITP reading (for consistency with other methods)
    void read_psf(const std::string& psf_file) { load_psf(psf_file); }
    void read_itp(const std::string& itp_file) { load_itp(itp_file); }

    // New methods from System class
    void load_structure_psf(const std::string& pdb, const std::string& psf);
    void load_structure_psf(const std::string& pdb, const std::vector<std::string>& psf_files);
    void load_structure_top(const std::string& pdb, const std::string& top);
    void load_structure_from_kwargs(const std::unordered_map<std::string, std::variant<std::string, std::vector<std::string>>>& kwargs);
    void load_structure_psf_auto(const std::string& pdb_file, const std::string& psf_file);
    void load_structure_psf_multi(const std::string& pdb_file, const std::string& psf_file);
    void load_structure_psf_single(const std::string& pdb_file, const std::string& psf_file, const std::string& target_residue);
    
    // PDB atom management methods from System
    size_t get_pdb_atom_count() const { return atoms_.size(); }
    const io::PDBAtom& get_pdb_atom(size_t index) const;
    io::PDBAtom& get_pdb_atom(size_t index);
    void add_pdb_atom(const io::PDBAtom& atom);
    void remove_pdb_atom(size_t index);
    std::vector<io::PDBAtom> get_pdb_atoms_by_residue(const std::string& residue_name) const;
    std::vector<io::PDBAtom> get_pdb_atoms_by_residue_sequence(const std::string& residue_name, int sequence) const;
    std::vector<io::PDBAtom> get_pdb_atoms_by_chain(char chain) const;
    void clear_pdb_atoms();
    bool has_pdb_atoms() const { return !atoms_.empty(); }

    // Static factory methods from System
    static Structure from_pdb_psf(const std::string& pdb_file, const std::vector<std::string>& psf_files);
    static Structure from_pdb_psf_itp(const std::string& pdb_file, const std::vector<std::string>& psf_files, const std::string& itp_file);
    static Structure from_pdb_psf_itps(const std::string& pdb_file, const std::vector<std::string>& psf_files, const std::vector<std::string>& itp_files);
    static Structure from_kwargs(const std::unordered_map<std::string, std::variant<std::string, std::vector<std::string>>>& kwargs);

private:
    std::vector<std::shared_ptr<io::IOResidue>> residues_;
    std::vector<std::shared_ptr<io::PDBAtom>> atoms_;
    std::optional<std::array<double, 6>> box_;  // Box dimensions (a, b, c, alpha, beta, gamma)
    std::shared_ptr<ForceField> forcefield_;
    
    // Cache for topology information
    std::optional<io::TopParser> cached_topology_;
    bool has_cached_topology_ = false;

    // Cache for PSF and ITP information
    std::optional<io::PSFParser> cached_psf_;
    std::vector<io::ITPParser> cached_itps_;
    bool has_cached_psf_ = false;

    // Helper functions for topology loading
    void update_atoms_topology(io::TopParser& top_parser);
    void update_atoms_topology(io::PSFParser& psf_parser);
    void update_atoms_topology(io::ITPParser& itp_parser);
    
    // Helper function to apply cached topology if exists
    void apply_cached_topology();
    void apply_cached_psf();
    void apply_cached_itps();

    // Private implementation functions
    void read_psf_file(const std::string& psf_file);
    void read_itp_file(const std::string& itp_file);
    
    // Helper function to validate indices
    void validate_pdb_atom_index(size_t index) const;
};

} // namespace core
} // namespace pygcmc 