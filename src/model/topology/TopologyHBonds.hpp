#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_HBONDS_HPP
#define PYGCMC_MODEL_TOPOLOGY_HBONDS_HPP

#include "TopologyCore.hpp"
#include "TopologyValidation.hpp"
#include <stdexcept>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Hydrogen bond management operations (donors and acceptors)
 */
class TopologyHBondManager : public ValidationMixin {
public:
    TopologyHBondManager(TopologyStorage& storage) : ValidationMixin(storage), storage_(storage) {}

    /**
     * @brief Add a hydrogen bond donor
     */
    void add_donor(int donor, int hydrogen) {
        if (!are_valid_bond_atoms(donor, hydrogen)) {
            throw std::invalid_argument("Invalid atom indices for donor");
        }
        
        TopologyDonor d;
        d.donor_atom = donor;
        d.hydrogen_atom = hydrogen;
        storage_.donors.push_back(d);
    }

    /**
     * @brief Add a hydrogen bond acceptor
     */
    void add_acceptor(int acceptor) {
        if (!is_valid_atom(acceptor)) {
            throw std::invalid_argument("Invalid atom index for acceptor");
        }
        
        TopologyAcceptor a;
        a.acceptor_atom = acceptor;
        storage_.acceptors.push_back(a);
    }

    // Getters
    const std::vector<TopologyDonor>& get_donors() const { return storage_.donors; }
    const std::vector<TopologyAcceptor>& get_acceptors() const { return storage_.acceptors; }

    // Count methods
    size_t get_num_donors() const { return storage_.donors.size(); }
    size_t get_num_acceptors() const { return storage_.acceptors.size(); }

    // Check methods
    bool has_donor(int donor_atom) const {
        return std::any_of(storage_.donors.begin(), storage_.donors.end(),
            [donor_atom](const TopologyDonor& donor) {
                return donor.donor_atom == donor_atom;
            });
    }

    bool has_donor(int donor_atom, int hydrogen_atom) const {
        return std::any_of(storage_.donors.begin(), storage_.donors.end(),
            [donor_atom, hydrogen_atom](const TopologyDonor& donor) {
                return donor.donor_atom == donor_atom && donor.hydrogen_atom == hydrogen_atom;
            });
    }

    bool has_acceptor(int acceptor_atom) const {
        return std::any_of(storage_.acceptors.begin(), storage_.acceptors.end(),
            [acceptor_atom](const TopologyAcceptor& acceptor) {
                return acceptor.acceptor_atom == acceptor_atom;
            });
    }

private:
    TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_HBONDS_HPP 