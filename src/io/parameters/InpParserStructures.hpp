#pragma once
#ifndef PYGCMC_IO_PARAMETERS_INPPARSERSTRUCTURES_HPP
#define PYGCMC_IO_PARAMETERS_INPPARSERSTRUCTURES_HPP

#include <string>
#include <vector>
#include <array>

namespace pygcmc {
namespace io {
namespace parameters {

/**
 * @brief Utility functions for INP file parsing
 */
class InpParserStructures {
public:
    /**
     * @brief Split string by delimiter
     * @param str Input string to split
     * @param delim Delimiter character (default: space)
     * @return std::vector<std::string> Vector of split tokens
     */
    static std::vector<std::string> split(const std::string& str, char delim = ' ');

    /**
     * @brief Trim whitespace from string
     * @param str Input string to trim
     * @return std::string Trimmed string
     */
    static std::string trim(const std::string& str);

    /**
     * @brief Parse array of 3 floats from string
     * @param str Input string containing 3 space-separated float values
     * @return std::array<float, 3> Array of parsed float values
     * @throw std::runtime_error If parsing fails
     */
    static std::array<float, 3> parse_float_array(const std::string& str);

    /**
     * @brief Parse vector of strings from string
     * @param str Input string containing space-separated string values
     * @return std::vector<std::string> Vector of parsed string values
     */
    static std::vector<std::string> parse_string_vector(const std::string& str);

    /**
     * @brief Parse vector of floats from string
     * @param str Input string containing space-separated float values
     * @return std::vector<float> Vector of parsed float values
     */
    static std::vector<float> parse_float_vector(const std::string& str);

    /**
     * @brief Parse vector of integers from string
     * @param str Input string containing space-separated integer values
     * @return std::vector<int> Vector of parsed integer values
     */
    static std::vector<int> parse_int_vector(const std::string& str);
};

} // namespace parameters
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_PARAMETERS_INPPARSERSTRUCTURES_HPP
