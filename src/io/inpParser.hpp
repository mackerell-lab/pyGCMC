#pragma once
#ifndef PYGCMC_IO_INPPARSER_HPP
#define PYGCMC_IO_INPPARSER_HPP

#include <string>
#include "model/param.hpp"

namespace pygcmc {
namespace io {

/**
 * @brief GCMC input file parser class
 * @details Used to parse GCMC simulation input files and store parameters in the Param data structure
 * Supported parameters include:
 * - File path parameters: par, fragitp, atomtypes, top, pdb, etc.
 * - Spatial parameters: grid_dx, box_size, cutoff, etc.
 * - Fragment parameters: fragname, fragconc, fragmuex
 * - Simulation control parameters: nprint, mcsteps, etc.
 * - Biased sampling parameters: use_cavity_bias, use_conf_bias
 */
class INPParser {
public:
    /**
     * @brief Parse input file and return a new Param object
     * @param filename Input file path
     * @return model::Param Param object containing parsing results
     * @throw std::runtime_error If the file doesn't exist or parsing fails
     */
    static model::Param parse_file(const std::string& filename);

    /**
     * @brief Parse input string and return a new Param object
     * @param content Input file content string
     * @return model::Param Param object containing parsing results
     * @throw std::runtime_error If parsing fails
     */
    static model::Param parse_string(const std::string& content);

    /**
     * @brief Parse input file and store results in an existing Param object
     * @param filename Input file path
     * @param param Param object for storing results
     * @throw std::runtime_error If the file doesn't exist or parsing fails
     */
    static void parse_to_param(const std::string& filename, model::Param& param);

    /**
     * @brief Parse input string and store results in an existing Param object
     * @param content Input file content string
     * @param param Param object for storing results
     * @throw std::runtime_error If parsing fails
     */
    static void parse_string_to_param(const std::string& content, model::Param& param);

private:
    /**
     * @brief Parse single line parameter
     * @param key Parameter name
     * @param value Parameter value
     * @param param Param object for storing results
     * @throw std::runtime_error If parsing fails
     */
    static void parse_line(const std::string& key, const std::string& value, model::Param& param);

    /**
     * @brief Validate parameter validity and consistency
     * @param param Param object to be validated
     * @throw std::runtime_error If validation fails
     */
    static void validate_parameters(model::Param& param);
};

} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_INPPARSER_HPP
