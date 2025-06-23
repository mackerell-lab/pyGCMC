#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_ATOMS_HPP
#define PYGCMC_MODEL_TOPOLOGY_ATOMS_HPP

#include "TopologyCore.hpp"
#include "TopologyValidationMixin.hpp"
#include <stdexcept>
#include <iostream>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Atom, residue, and segment management operations
 */
class TopologyAtomManager : public ValidationMixin {
public:
    TopologyAtomManager(TopologyStorage& storage) : ValidationMixin(storage), storage_(storage) {}

    /**
     * @brief Add an atom to the topology
     * @param name Atom name (e.g., "CA", "N", "O")
     * @param type Atom type (e.g., "CT1", "NH1")
     * @param charge Atomic charge
     * @param mass Atomic mass
     * @param residue_name Residue name
     * @param residue_number Residue number
     * @param segment_name Segment name
     * @return Atom ID
     */
    int add_atom(const std::string& name, const std::string& type, double charge, double mass,
                const std::string& residue_name, int residue_number, const std::string& segment_name) {
        try {
            // First ensure we have the segment
            int segment_id = ensure_segment(segment_name);
            
            // Then ensure we have the residue
            int residue_id = ensure_residue(residue_name, residue_number, segment_name, segment_id);

            // Create and add the atom
            TopologyAtom atom;
            atom.id = static_cast<int>(storage_.atoms.size());
            atom.name = name;
            atom.type = type;
            atom.charge = charge;
            atom.mass = mass;
            atom.residue_id = residue_id;
            atom.segment_id = segment_id;

            // Add atom to the vectors and maps
            int atom_id = static_cast<int>(storage_.atoms.size());
            storage_.atoms.push_back(atom);
            storage_.atom_map[std::make_tuple(residue_name, residue_number, segment_name, name)] = atom_id;

            // Add atom to its residue
            if (is_valid_residue(residue_id)) {
                storage_.residues[residue_id].atoms.push_back(atom_id);
            } else {
                throw std::runtime_error("Invalid residue ID");
            }

            return atom_id;
        } catch (const std::exception& e) {
            std::cerr << "Error in add_atom: " << e.what() << std::endl;
            throw;
        }
    }

    /**
     * @brief Add a residue to the topology
     * @param name Residue name
     * @param number Residue number
     * @param segment Segment name
     * @return Residue ID
     */
    int add_residue(const std::string& name, int number, const std::string& segment) {
        // First ensure we have the segment
        int segment_id = ensure_segment(segment);

        // Create and add the residue
        TopologyResidue residue;
        residue.id = static_cast<int>(storage_.residues.size());
        residue.name = name;
        residue.number = number;
        residue.segment = segment;

        // Add residue to the vectors and maps
        int residue_id = static_cast<int>(storage_.residues.size());
        storage_.residues.push_back(residue);
        storage_.residue_map[std::make_tuple(name, number, segment)] = residue_id;

        // Add residue to its segment
        if (is_valid_segment(segment_id)) {
            storage_.segments[segment_id].residues.push_back(residue_id);
        }

        return residue_id;
    }

    /**
     * @brief Add a segment to the topology
     * @param name Segment name
     * @return Segment ID
     */
    int add_segment(const std::string& name) {
        TopologySegment segment;
        segment.id = static_cast<int>(storage_.segments.size());
        segment.name = name;

        // Add segment to the vectors and maps
        int segment_id = static_cast<int>(storage_.segments.size());
        storage_.segments.push_back(segment);
        storage_.segment_map[name] = segment_id;

        return segment_id;
    }

    /**
     * @brief Get atom by index
     */
    const TopologyAtom& get_atom(int index) const {
        if (!is_valid_atom(index)) {
            throw std::out_of_range("Invalid atom index");
        }
        return storage_.atoms[index];
    }

    /**
     * @brief Get residue by index
     */
    const TopologyResidue& get_residue(int index) const {
        if (!is_valid_residue(index)) {
            throw std::out_of_range("Invalid residue index");
        }
        return storage_.residues[index];
    }

    /**
     * @brief Get segment by index
     */
    const TopologySegment& get_segment(int index) const {
        if (!is_valid_segment(index)) {
            throw std::out_of_range("Invalid segment index");
        }
        return storage_.segments[index];
    }

    /**
     * @brief Get number of atoms
     */
    int get_num_atoms() const { 
        return static_cast<int>(storage_.atoms.size()); 
    }

    /**
     * @brief Get number of residues
     */
    int get_num_residues() const { 
        return static_cast<int>(storage_.residues.size()); 
    }

    /**
     * @brief Get number of segments
     */
    int get_num_segments() const { 
        return static_cast<int>(storage_.segments.size()); 
    }

    /**
     * @brief Find atom by residue and atom name
     */
    std::optional<int> find_atom(const std::string& residue_name, int residue_number,
                                const std::string& atom_name) const {
        // Try to find the atom in any segment
        for (const auto& segment : storage_.segments) {
            auto it = storage_.atom_map.find(std::make_tuple(residue_name, residue_number, segment.name, atom_name));
            if (it != storage_.atom_map.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Find residue by name and number
     */
    std::optional<int> find_residue(const std::string& name, int number) const {
        // Try to find the residue in any segment
        for (const auto& segment : storage_.segments) {
            auto it = storage_.residue_map.find(std::make_tuple(name, number, segment.name));
            if (it != storage_.residue_map.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Find segment by name
     */
    std::optional<int> find_segment(const std::string& name) const {
        auto it = storage_.segment_map.find(name);
        if (it != storage_.segment_map.end()) {
            return it->second;
        }
        return std::nullopt;
    }

    /**
     * @brief Reserve space for atoms
     */
    void reserve_atoms(size_t n) {
        storage_.atoms.reserve(n);
    }

private:
    TopologyStorage& storage_;

    /**
     * @brief Ensure segment exists, create if necessary
     */
    int ensure_segment(const std::string& segment_name) {
        auto segment_it = storage_.segment_map.find(segment_name);
        if (segment_it == storage_.segment_map.end()) {
            return add_segment(segment_name);
        }
        return segment_it->second;
    }

    /**
     * @brief Ensure residue exists, create if necessary
     */
    int ensure_residue(const std::string& residue_name, int residue_number, 
                      const std::string& segment_name, int segment_id) {
        auto residue_key = std::make_tuple(residue_name, residue_number, segment_name);
        auto residue_it = storage_.residue_map.find(residue_key);
        if (residue_it == storage_.residue_map.end()) {
            // Create new residue
            TopologyResidue residue;
            residue.id = static_cast<int>(storage_.residues.size());
            residue.name = residue_name;
            residue.number = residue_number;
            residue.segment = segment_name;
            int residue_id = residue.id;
            storage_.residues.push_back(residue);
            storage_.residue_map[residue_key] = residue_id;
            
            // Add residue to its segment
            if (is_valid_segment(segment_id)) {
                storage_.segments[segment_id].residues.push_back(residue_id);
            } else {
                throw std::runtime_error("Invalid segment ID");
            }
            
            return residue_id;
        }
        return residue_it->second;
    }
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_ATOMS_HPP 