// src/model/forcefield.hpp

#pragma once

// Include the new refactored forcefield module
#include "topology/ForceFieldMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using ForceField = topology::ForceField;
using NonbondedParams = topology::NonbondedParams;
using LJParams = topology::LJParams;
using BondParams = topology::BondParams;
using AngleParams = topology::AngleParams;
using DihedralParams = topology::DihedralParams;
using ImproperParams = topology::ImproperParams;
using NBFIXParams = topology::NBFIXParams;

} // namespace model
} // namespace pygcmc


