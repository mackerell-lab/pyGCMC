// src/model/structure.hpp

#pragma once
#ifndef PYGCMC_MODEL_STRUCTURE_HPP
#define PYGCMC_MODEL_STRUCTURE_HPP

// Include the new refactored structure module
#include "structure/StructureMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using Structure = structure::Structure;

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_STRUCTURE_HPP
