#pragma once

#ifndef PYGCMC_MODEL_PARAM_QUERIES_HPP
#define PYGCMC_MODEL_PARAM_QUERIES_HPP

#include "ParamStructures.hpp"
#include <sstream>
#include <stdexcept>

namespace pygcmc {
namespace model {
namespace param {

/**
 * @brief Static query and validation operations for parameters
 */
class ParamQueries {
public:
    // === Validation methods ===
    static bool isValidMCParams(const MCParams& mc_params) {
        return mc_params.temperature > 0.0f && mc_params.mc_steps > 0;
    }

    static bool isValidEnergyInfo(const EnergyInfo& energy_info) {
        return energy_info.fragment_cutoff > 0.0f && energy_info.protein_cutoff > 0.0f;
    }

    static bool hasRequiredFiles(const FileInfo& file_info) {
        return !file_info.topology_file.empty() && 
               !file_info.input_pdb_file.empty() && 
               !file_info.output_pdb_file.empty();
    }

    static bool isValidFileInfo(const FileInfo& file_info) {
        return hasRequiredFiles(file_info);
    }

    static bool isValidParam(const BasicInfo& basic_info, const SpaceInfo& space_info, 
                           const MCParams& mc_params, const EnergyInfo& energy_info,
                           const FragmentInfo& fragment_info, const BiasInfo& bias_info,
                           const FileInfo& file_info) {
        (void)basic_info;
        (void)space_info;
        (void)fragment_info;
        (void)bias_info;
        return isValidMCParams(mc_params) && 
               isValidEnergyInfo(energy_info) && 
               isValidFileInfo(file_info);
    }

    // === Query helper methods ===
    static std::string getFragmentCoordinateFilename(const FileInfo& file_info, int frag_index) {
        if (frag_index < 0 || frag_index >= static_cast<int>(file_info.fragment_names.size())) {
            throw std::out_of_range("Fragment index out of range");
        }
        return file_info.monomer_dir + "/" + file_info.fragment_names[frag_index] + ".pdb";
    }

    // === String representation ===
    static std::string toString(const BasicInfo& basic_info, const SpaceInfo& space_info,
                              const MCParams& mc_params, const EnergyInfo& energy_info,
                              const FragmentInfo& fragment_info, const BiasInfo& bias_info,
                              const FileInfo& file_info) {
        (void)energy_info;
        (void)fragment_info;
        (void)bias_info;
        std::stringstream ss;
        ss << "GCMC Parameters:\n";
        ss << "  Version: " << basic_info.version << "\n";
        ss << "  Temperature: " << mc_params.temperature << " K\n";
        ss << "  MC Steps: " << mc_params.mc_steps << "\n";
        ss << "  Box Size: [" << space_info.box_size[0] << ", " 
           << space_info.box_size[1] << ", " << space_info.box_size[2] << "] Å\n";
        ss << "  Cutoff: " << space_info.cutoff << " Å\n";
        ss << "  Number of fragments: " << file_info.fragment_names.size();
        return ss.str();
    }

    // === Component-specific queries ===
    static std::string getBasicInfoString(const BasicInfo& basic_info) {
        std::stringstream ss;
        ss << "Basic Info: version=" << basic_info.version 
           << ", verbosity=" << basic_info.verbosity
           << ", debug=" << basic_info.debug;
        return ss.str();
    }

    static std::string getSpaceInfoString(const SpaceInfo& space_info) {
        std::stringstream ss;
        ss << "Space Info: box=[" << space_info.box_size[0] << "," 
           << space_info.box_size[1] << "," << space_info.box_size[2] 
           << "], volume=" << space_info.volume << ", cutoff=" << space_info.cutoff;
        return ss.str();
    }

    static std::string getMCParamsString(const MCParams& mc_params) {
        std::stringstream ss;
        ss << "MC Params: T=" << mc_params.temperature << "K, steps=" << mc_params.mc_steps
           << ", beta=" << mc_params.beta;
        return ss.str();
    }
};

} // namespace param
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_PARAM_QUERIES_HPP