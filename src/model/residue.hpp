// src/model/residue.hpp

#pragma once
#ifndef PYGCMC_MODEL_RESIDUE_HPP
#define PYGCMC_MODEL_RESIDUE_HPP

// Include the new refactored residue module
#include "residue/ResidueMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using Residue = residue::Residue;

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_RESIDUE_HPP