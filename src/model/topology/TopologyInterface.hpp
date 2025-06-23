#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_INTERFACE_HPP
#define PYGCMC_MODEL_TOPOLOGY_INTERFACE_HPP

#include "../common/ModelInterface.hpp"
#include "TopologyCore.hpp"
#include <vector>
#include <string>
#include <optional>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Interface for topology data access
 */
class ITopologyReader {
public:
    virtual ~ITopologyReader() = default;

    // Basic getters
    virtual const TopologyAtom& get_atom(int index) const = 0;
    virtual const TopologyResidue& get_residue(int index) const = 0;
    virtual const TopologySegment& get_segment(int index) const = 0;
    virtual const std::vector<TopologyBond>& get_bonds() const = 0;
    virtual const std::vector<TopologyAngle>& get_angles() const = 0;
    virtual const std::vector<TopologyDihedral>& get_dihedrals() const = 0;

    // Counts
    virtual int get_num_atoms() const = 0;
    virtual int get_num_residues() const = 0;
    virtual int get_num_segments() const = 0;
    virtual size_t get_num_bonds() const = 0;
    virtual size_t get_num_angles() const = 0;
    virtual size_t get_num_dihedrals() const = 0;
    virtual size_t get_num_impropers() const = 0;

    // Existence checks
    virtual bool has_atom(int index) const = 0;
    virtual bool has_residue(int index) const = 0;
    virtual bool has_segment(int index) const = 0;

    // Search functions
    virtual std::optional<int> find_atom(const std::string& residue_name, int residue_number,
                                        const std::string& atom_name) const = 0;
    virtual std::optional<int> find_residue(const std::string& name, int number) const = 0;
    virtual std::optional<int> find_segment(const std::string& name) const = 0;
};

/**
 * @brief Interface for topology modification
 */
class ITopologyWriter {
public:
    virtual ~ITopologyWriter() = default;

    // Add structural elements
    virtual int add_atom(const std::string& name, const std::string& type, double charge, double mass,
                        const std::string& residue_name, int residue_number, const std::string& segment_name) = 0;
    virtual int add_residue(const std::string& name, int number, const std::string& segment) = 0;
    virtual int add_segment(const std::string& name) = 0;

    // Add connectivity
    virtual void add_bond(int atom1, int atom2, double length = 0.0, double force_constant = 0.0, int function_type = 1) = 0;
    virtual void add_angle(int atom1, int atom2, int atom3, double angle = 0.0, double force_constant = 0.0, int function_type = 1) = 0;
    virtual void add_dihedral(int atom1, int atom2, int atom3, int atom4, int multiplicity = 1,
                             double angle = 0.0, double force_constant = 0.0, bool improper = false, int function_type = 1) = 0;

    // Utility
    virtual void reserve_atoms(size_t n) = 0;
};

/**
 * @brief Interface for advanced topology features
 */
class ITopologyAdvanced {
public:
    virtual ~ITopologyAdvanced() = default;

    // Advanced features
    virtual void add_title(const std::string& title) = 0;
    virtual void add_improper(int atom1, int atom2, int atom3, int atom4, double angle = 0.0, double force_constant = 0.0) = 0;
    virtual void add_donor(int donor, int hydrogen) = 0;
    virtual void add_acceptor(int acceptor) = 0;
    virtual void add_nonbonded_exclusion(int atom1, int atom2) = 0;
    virtual void add_group(int id, const std::vector<int>& atoms, const std::string& type = "") = 0;
    virtual void add_cmap(const std::array<int, 8>& atoms) = 0;
    virtual void add_cmap(const std::array<int, 5>& atoms, int function_type = 1) = 0;

    // Advanced getters
    virtual const std::vector<std::string>& get_titles() const = 0;
    virtual const std::vector<TopologyDonor>& get_donors() const = 0;
    virtual const std::vector<TopologyAcceptor>& get_acceptors() const = 0;
    virtual const std::vector<TopologyGroup>& get_groups() const = 0;
    virtual const std::vector<TopologyCmap>& get_cmaps() const = 0;
    virtual const std::map<int, std::set<int>>& get_exclusions() const = 0;

    // Advanced checks
    virtual bool has_donor(int donor_atom) const = 0;
    virtual bool has_donor(int donor_atom, int hydrogen_atom) const = 0;
    virtual bool has_acceptor(int acceptor_atom) const = 0;
    virtual bool has_cmap() const = 0;
    virtual bool has_cmap(const std::vector<int>& atoms) const = 0;
    virtual bool has_group(int group_id) const = 0;
    virtual const TopologyGroup& get_group(int index) const = 0;

    // Connectivity checks
    virtual bool has_bond(int atom1, int atom2) const = 0;
    virtual bool has_angle(int atom1, int atom2, int atom3) const = 0;
    virtual bool has_dihedral(int atom1, int atom2, int atom3, int atom4) const = 0;
    virtual bool has_improper(int atom1, int atom2, int atom3, int atom4) const = 0;

    // Counts for advanced features
    virtual size_t get_num_donors() const = 0;
    virtual size_t get_num_acceptors() const = 0;
    virtual size_t get_num_cmaps() const = 0;
    virtual size_t get_num_groups() const = 0;
};

/**
 * @brief Complete topology interface combining all capabilities
 */
class ITopology : public ITopologyReader, 
                  public ITopologyWriter, 
                  public ITopologyAdvanced,
                  public common::IValidatable {
public:
    virtual ~ITopology() = default;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_INTERFACE_HPP 