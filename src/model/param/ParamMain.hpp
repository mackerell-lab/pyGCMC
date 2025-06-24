#pragma once

#ifndef PYGCMC_MODEL_PARAM_MAIN_HPP
#define PYGCMC_MODEL_PARAM_MAIN_HPP

#include "../common/ModelInterface.hpp"
#include "ParamStructures.hpp"
#include "ParamOperations.hpp"
#include "ParamQueries.hpp"
#include <memory>

namespace pygcmc {
namespace model {
namespace param {

/**
 * @brief Complete GCMC parameter class with modular design and backward compatibility
 * Uses delegation pattern to maintain clean separation of concerns
 */
class Param : public common::IValidatable {
public:
    Param() = default;
    ~Param() = default;

    // === Backward Compatibility Type Aliases ===
    using BasicInfo = ::pygcmc::model::param::BasicInfo;
    using SpaceInfo = ::pygcmc::model::param::SpaceInfo;
    using MCInfo = ::pygcmc::model::param::MCParams;
    using EnergyInfo = ::pygcmc::model::param::EnergyInfo;
    using FragmentInfo = ::pygcmc::model::param::FragmentInfo;
    using BiasInfo = ::pygcmc::model::param::BiasInfo;
    using FileInfo = ::pygcmc::model::param::FileInfo;

    // === IValidatable Interface ===
    bool is_valid() const override {
        return ParamQueries::isValidParam(basic_info_, space_info_, mc_info_, 
                                        energy_info_, fragment_info_, bias_info_, file_info_);
    }

    // === Const Getters ===
    const BasicInfo& get_basic_info() const { return basic_info_; }
    const SpaceInfo& get_space_info() const { return space_info_; }
    const MCInfo& get_mc_info() const { return mc_info_; }
    const EnergyInfo& get_energy_info() const { return energy_info_; }
    const FragmentInfo& get_fragment_info() const { return fragment_info_; }
    const BiasInfo& get_bias_info() const { return bias_info_; }
    const FileInfo& get_file_info() const { return file_info_; }

    // === Non-const Getters ===
    BasicInfo& get_basic_info() { return basic_info_; }
    SpaceInfo& get_space_info() { return space_info_; }
    MCInfo& get_mc_info() { return mc_info_; }
    EnergyInfo& get_energy_info() { return energy_info_; }
    FragmentInfo& get_fragment_info() { return fragment_info_; }
    BiasInfo& get_bias_info() { return bias_info_; }
    FileInfo& get_file_info() { return file_info_; }

    // === Setters ===
    void set_basic_info(const BasicInfo& info) { basic_info_ = info; }
    void set_space_info(const SpaceInfo& info) { space_info_ = info; }
    void set_mc_info(const MCInfo& info) { mc_info_ = info; }
    void set_energy_info(const EnergyInfo& info) { energy_info_ = info; }
    void set_fragment_info(const FragmentInfo& info) { fragment_info_ = info; }
    void set_bias_info(const BiasInfo& info) { bias_info_ = info; }
    void set_file_info(const FileInfo& info) { file_info_ = info; }

    // === Convenience Methods (Delegate to Operations) ===
    void set_temperature(float temperature) {
        ParamOperations::setTemperature(mc_info_, temperature);
    }

    void set_box_size(float x, float y, float z) {
        ParamOperations::setBoxSize(space_info_, x, y, z);
    }

    void set_mc_steps(int steps) {
        ParamOperations::setMCSteps(mc_info_, steps);
    }

    void set_cutoff(float cutoff) {
        ParamOperations::setCutoff(space_info_, energy_info_, cutoff);
    }

    // === Update Operations (Delegate to Operations) ===
    void update_derived_values() {
        ParamOperations::updateBeta(mc_info_);
        ParamOperations::updateVolume(space_info_);
        ParamOperations::updateEnergySquaredValues(energy_info_);
        ParamOperations::updateFragmentSquaredValues(fragment_info_);
        ParamOperations::updateBiasSquaredValues(bias_info_);
    }

    // === Clear Operations (Delegate to Operations) ===
    void clear() {
        ParamOperations::clearBasicInfo(basic_info_);
        ParamOperations::clearSpaceInfo(space_info_);
        ParamOperations::clearMCParams(mc_info_);
        ParamOperations::clearEnergyInfo(energy_info_);
        ParamOperations::clearFragmentInfo(fragment_info_);
        ParamOperations::clearBiasInfo(bias_info_);
        ParamOperations::clearFileInfo(file_info_);
    }

    // === Query Methods (Delegate to Queries) ===
    std::string get_fragment_coordinate_filename(int frag_index) const {
        return ParamQueries::getFragmentCoordinateFilename(file_info_, frag_index);
    }

    std::string to_string() const {
        return ParamQueries::toString(basic_info_, space_info_, mc_info_, 
                                    energy_info_, fragment_info_, bias_info_, file_info_);
    }

    // === Utility Methods ===
    std::unique_ptr<Param> clone() const {
        auto cloned = std::make_unique<Param>();
        cloned->basic_info_ = basic_info_;
        cloned->space_info_ = space_info_;
        cloned->mc_info_ = mc_info_;
        cloned->energy_info_ = energy_info_;
        cloned->fragment_info_ = fragment_info_;
        cloned->bias_info_ = bias_info_;
        cloned->file_info_ = file_info_;
        return cloned;
    }

private:
    BasicInfo basic_info_;
    SpaceInfo space_info_;
    MCInfo mc_info_;
    EnergyInfo energy_info_;
    FragmentInfo fragment_info_;
    BiasInfo bias_info_;
    FileInfo file_info_;
};

} // namespace param
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_MAIN_HPP