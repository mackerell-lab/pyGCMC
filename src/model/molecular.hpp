#pragma once
#ifndef PYGCMC_MODEL_MOLECULAR_HPP
#define PYGCMC_MODEL_MOLECULAR_HPP

// Include the new refactored molecular module
#include "molecule/MolecularMain.hpp"
// Include necessary dependencies for type aliases
#include "atom.hpp"
#include "residue.hpp"
#include "topology.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using Molecular = molecule::Molecular;

// Re-export commonly used types for system compatibility
using Atom = atom::Atom;
using Residue = residue::Residue;
using Topology = topology::Topology;

// Re-export topology types that might be used by system files
using TopologyAtom = topology::TopologyAtom;
using TopologyResidue = topology::TopologyResidue;
using TopologyBond = topology::TopologyBond;
using TopologyAngle = topology::TopologyAngle;
using TopologyDihedral = topology::TopologyDihedral;
using TopologyCmap = topology::TopologyCmap;
using TopologyDonor = topology::TopologyDonor;
using TopologyAcceptor = topology::TopologyAcceptor;
using TopologyGroup = topology::TopologyGroup;
using TopologySegment = topology::TopologySegment;

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULAR_HPP
