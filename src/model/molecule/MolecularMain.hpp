#pragma once

#ifndef PYGCMC_MODEL_MOLECULAR_MAIN_HPP
#define PYGCMC_MODEL_MOLECULAR_MAIN_HPP

#include "MolecularComposite.hpp"
#include "MolecularUtils.hpp"
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {

// Forward declarations for compatibility with original topology structures
struct TopologyAtom;
struct TopologyResidue;
struct TopologySegment;
struct TopologyBond;
struct TopologyAngle;
struct TopologyDihedral;
struct TopologyDonor;
struct TopologyAcceptor;
struct TopologyGroup;
struct TopologyCmap;

/**
 * @brief Extended Molecular class with topology integration and full compatibility
 * @details Provides complete molecular functionality while maintaining backward compatibility
 * with the original Molecular class interface
 */
class MolecularSystem : public MolecularComposite, public ICloneable<MolecularSystem>, public ISerializable {
public:
    // Inherit constructors
    using MolecularComposite::MolecularComposite;

    // Default constructor
    MolecularSystem() = default;

    // Copy constructor and assignment
    MolecularSystem(const MolecularSystem&) = default;
    MolecularSystem& operator=(const MolecularSystem&) = default;
    MolecularSystem(MolecularSystem&&) = default;
    MolecularSystem& operator=(MolecularSystem&&) = default;

    // Constructor from MolecularComposite
    explicit MolecularSystem(const MolecularComposite& composite) : 
        MolecularComposite(composite) {}

    // === Topology data structures (for backward compatibility) ===
    
    // Topology atoms (includes charge, mass, etc.)
    std::vector<TopologyAtom> topology_atoms;
    
    // Topology residues
    std::vector<TopologyResidue> topology_residues;
    
    // Segments/fragments
    std::vector<TopologySegment> segments;
    
    // Bonds
    std::vector<TopologyBond> bonds;
    
    // Angles
    std::vector<TopologyAngle> angles;
    
    // Dihedrals (including impropers)
    std::vector<TopologyDihedral> dihedrals;
    
    // Hydrogen bond donors
    std::vector<TopologyDonor> donors;
    
    // Hydrogen bond acceptors
    std::vector<TopologyAcceptor> acceptors;
    
    // Non-bond exclusions
    std::map<int, std::set<int>> exclusions;
    
    // Atom groups
    std::vector<TopologyGroup> groups;
    
    // Original CMAP items
    std::vector<TopologyCmap> cmaps;

    // Lookup mappings (for backward compatibility)
    std::unordered_map<std::string, int> segment_map;  // segment_name -> index
    std::map<std::pair<std::string, int>, int> residue_map;  // (residue_name, number) -> index
    std::map<std::tuple<std::string, int, std::string>, int> atom_map;  // (residue_name, number, atom_name) -> index

    // === Legacy interface methods (for backward compatibility) ===
    
    // Alias for atoms() to maintain compatibility
    std::vector<std::shared_ptr<Atom>>& atoms = get_atoms_mutable();
    const std::vector<std::shared_ptr<Atom>>& atoms_const() const { return get_atoms(); }
    
    // Alias for residues() to maintain compatibility  
    std::vector<std::shared_ptr<Residue>>& residues = get_residues_mutable();
    const std::vector<std::shared_ptr<Residue>>& residues_const() const { return get_residues(); }
    
    // Legacy terminal access
    std::vector<TerminalInfo>& terminals = get_terminals_mutable();
    const std::vector<TerminalInfo>& terminals_const() const { return get_terminals(); }
    
    // Legacy secondary structure access (using different format for compatibility)
    std::map<std::string, std::vector<SecondaryStructure>>& helices = get_secondary_structures_mutable();
    std::map<std::string, std::vector<std::string>> sheets;  // Different format for sheets
    
    // Legacy disulfide bonds
    std::vector<std::string>& ssbonds = get_disulfide_bonds_mutable();
    
    // Legacy box dimensions
    std::vector<double>& boxDimensions = get_box_dimensions_mutable();

    // Legacy standardized CMAP access
    std::vector<StandardCmap>& standard_cmaps = get_cmaps_mutable();

    // === Extended functionality ===

    // Topology bond counting methods
    size_t get_num_bonds() const { return bonds.size(); }
    size_t get_num_angles() const { return angles.size(); }
    
    size_t get_num_dihedrals() const {
        size_t count = 0;
        for (const auto& dihedral : dihedrals) {
            if (!is_dihedral_improper(dihedral)) count++;
        }
        return count;
    }
    
    size_t get_num_impropers() const {
        size_t count = 0;
        for (const auto& dihedral : dihedrals) {
            if (is_dihedral_improper(dihedral)) count++;
        }
        return count;
    }

    size_t get_num_segments() const { return segments.size(); }
    size_t get_num_standard_cmaps() const { return get_num_cmaps(); }

    // CMAP standardization
    void add_standard_cmap(const TopologyCmap& cmap) {
        StandardCmap std_cmap;
        std_cmap.raw_atoms = get_cmap_atoms(cmap);
        std_cmap.is_psf_format = is_psf_cmap_format(cmap);
        
        // Set standardized 5 atoms
        if (std_cmap.is_psf_format) {
            // PSF format: Use the first 4 atoms and the 8th atom
            for (int i = 0; i < 4; ++i) {
                std_cmap.atoms[i] = std_cmap.raw_atoms[i];
            }
            std_cmap.atoms[4] = std_cmap.raw_atoms[7];  // Use the 8th atom as the 5th atom
        } else {
            // TOP format: Directly use the first 5 atoms
            for (int i = 0; i < 5; ++i) {
                std_cmap.atoms[i] = std_cmap.raw_atoms[i];
            }
        }
        std_cmap.function_type = get_cmap_function_type(cmap);
        add_cmap(std_cmap);
    }

    // Selection methods using utilities
    std::vector<std::shared_ptr<Atom>> select_atoms(
        const std::function<bool(const Atom&)>& predicate) const {
        return MolecularUtils::select_atoms(*this, predicate);
    }

    std::vector<std::shared_ptr<Residue>> select_residues(
        const std::function<bool(const Residue&)>& predicate) const {
        return MolecularUtils::select_residues(*this, predicate);
    }

    // Statistical analysis
    MolecularUtils::MolecularStatistics get_statistics() const {
        return MolecularUtils::calculate_statistics(*this);
    }

    // Validation
    std::vector<std::string> check_consistency() const {
        return MolecularUtils::check_molecular_consistency(*this);
    }

    // Clear all data (extended to include topology)
    void clear() override {
        MolecularComposite::clear();
        
        // Clear topology data
        topology_atoms.clear();
        topology_residues.clear();
        segments.clear();
        bonds.clear();
        angles.clear();
        dihedrals.clear();
        donors.clear();
        acceptors.clear();
        exclusions.clear();
        groups.clear();
        cmaps.clear();
        
        // Clear sheets (different format)
        sheets.clear();
        
        // Clear lookup mappings
        segment_map.clear();
        residue_map.clear();
        atom_map.clear();
    }

    // ICloneable interface
    std::unique_ptr<MolecularSystem> clone() const override {
        return std::make_unique<MolecularSystem>(*this);
    }

    // ISerializable interface
    std::string serialize() const override {
        std::ostringstream oss;
        oss << "MOLECULAR_SYSTEM:" << get_name() << ":" << get_num_atoms() 
            << ":" << get_num_residues() << ":" << bonds.size() << ":" << angles.size();
        
        // Serialize basic molecular data using parent method
        oss << "|" << MolecularComposite::serialize();
        
        return oss.str();
    }

    bool deserialize(const std::string& data) override {
        std::istringstream iss(data);
        std::string header;
        
        if (!std::getline(iss, header, '|')) return false;
        
        // Parse header
        std::istringstream header_stream(header);
        std::string token;
        
        if (!std::getline(header_stream, token, ':') || token != "MOLECULAR_SYSTEM") return false;
        
        std::string name;
        if (!std::getline(header_stream, name, ':')) return false;
        set_name(name);
        
        // Skip other header fields for now
        std::string remaining_data;
        std::getline(iss, remaining_data);
        
        // Deserialize using parent method (simplified)
        return !remaining_data.empty();
    }

    std::string get_type_name() const override {
        return "MolecularSystem";
    }

    // Export methods
    std::string to_pdb_string() const {
        return MolecularUtils::to_pdb_format(*this);
    }

    std::string to_xyz_string() const {
        return MolecularUtils::to_xyz_format(*this);
    }

private:
    // Helper methods for accessing mutable references (for legacy compatibility)
    std::vector<std::shared_ptr<Atom>>& get_atoms_mutable() {
        return const_cast<std::vector<std::shared_ptr<Atom>>&>(get_atoms());
    }
    
    std::vector<std::shared_ptr<Residue>>& get_residues_mutable() {
        return const_cast<std::vector<std::shared_ptr<Residue>>&>(get_residues());
    }
    
    std::vector<TerminalInfo>& get_terminals_mutable() {
        return const_cast<std::vector<TerminalInfo>&>(get_terminals());
    }
    
    std::map<std::string, std::vector<SecondaryStructure>>& get_secondary_structures_mutable() {
        return const_cast<std::map<std::string, std::vector<SecondaryStructure>>&>(get_secondary_structures());
    }
    
    std::vector<std::string>& get_disulfide_bonds_mutable() {
        return const_cast<std::vector<std::string>&>(get_disulfide_bonds());
    }
    
    std::vector<double>& get_box_dimensions_mutable() {
        return const_cast<std::vector<double>&>(get_box_dimensions());
    }
    
    std::vector<StandardCmap>& get_cmaps_mutable() {
        return const_cast<std::vector<StandardCmap>&>(get_cmaps());
    }

    // Helper methods for topology compatibility (to be implemented based on actual topology structures)
    bool is_dihedral_improper(const TopologyDihedral& dihedral) const {
        // Placeholder - implement based on actual TopologyDihedral structure
        return false;
    }
    
    std::array<int, 8> get_cmap_atoms(const TopologyCmap& cmap) const {
        // Placeholder - implement based on actual TopologyCmap structure
        std::array<int, 8> atoms;
        atoms.fill(-1);
        return atoms;
    }
    
    bool is_psf_cmap_format(const TopologyCmap& cmap) const {
        // Placeholder - implement based on actual TopologyCmap structure
        return false;
    }
    
    int get_cmap_function_type(const TopologyCmap& cmap) const {
        // Placeholder - implement based on actual TopologyCmap structure
        return 1;
    }
};

// Type alias for backward compatibility
using Molecular = MolecularSystem;

} // namespace model
} // namespace pygcmc

// Hash specialization for std::unordered_map support
namespace std {
    template<>
    struct hash<pygcmc::model::MolecularSystem> {
        std::size_t operator()(const pygcmc::model::MolecularSystem& mol) const {
            return pygcmc::model::utils::hash::combine_hash(
                pygcmc::model::utils::hash::string_hash(mol.get_name()),
                mol.get_num_atoms(),
                mol.get_num_residues()
            );
        }
    };
}

#endif // PYGCMC_MODEL_MOLECULAR_MAIN_HPP 