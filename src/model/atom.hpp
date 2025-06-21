// src/model/atom.hpp

#pragma once

#ifndef PYGCMC_MODEL_ATOM_HPP
#define PYGCMC_MODEL_ATOM_HPP

// Include the new refactored atom module
#include "atom/AtomMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using Atom = atom::Atom;

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_ATOM_HPP

