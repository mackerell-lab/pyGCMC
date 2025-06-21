// src/model/param.hpp

#pragma once
#ifndef PYGCMC_MODEL_PARAM_HPP
#define PYGCMC_MODEL_PARAM_HPP

// Include the new refactored param module
#include "param/ParamMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using Param = param::Param;
using BasicInfo = param::BasicInfo;
using SpaceInfo = param::SpaceInfo;
using MCInfo = param::MCInfo;
using EnergyInfo = param::EnergyInfo;
using FragmentInfo = param::FragmentInfo;
using BiasInfo = param::BiasInfo;
using FileInfo = param::FileInfo;

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_HPP

 