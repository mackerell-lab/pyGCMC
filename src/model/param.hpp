// src/model/param.hpp

#pragma once
#ifndef PYGCMC_MODEL_PARAM_HPP
#define PYGCMC_MODEL_PARAM_HPP

// Include the new refactored param module
#include "param/ParamMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases - as standalone types
using BasicInfo = param::BasicInfo;
using SpaceInfo = param::SpaceInfo;
using EnergyInfo = param::EnergyInfo;
using FragmentInfo = param::FragmentInfo;
using BiasInfo = param::BiasInfo;
using FileInfo = param::FileInfo;

// Main Param class with nested type compatibility
class Param : public param::Param {
public:
    // Re-export types as nested types for Python binding compatibility
    using BasicInfo = param::BasicInfo;
    using SpaceInfo = param::SpaceInfo;
    using MCInfo = param::MCParams;  // Nested alias for backward compatibility
    using EnergyInfo = param::EnergyInfo;
    using FragmentInfo = param::FragmentInfo;
    using BiasInfo = param::BiasInfo;
    using FileInfo = param::FileInfo;
    
    // Inherit all constructors and methods
    using param::Param::Param;
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_HPP

 