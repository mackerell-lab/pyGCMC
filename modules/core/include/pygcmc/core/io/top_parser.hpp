// modules/core/include/pygcmc/core/io/top_parser.hpp

#ifndef PYGCMC_CORE_IO_TOP_PARSER_HPP
#define PYGCMC_CORE_IO_TOP_PARSER_HPP

#include "parser_common.hpp"
#include <string>
#include <istream>

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief TOP Parser
 * 
 * Parse TOP files to extract atom type information and force field parameters.
 */
class TopParser {
public:
    /**
     * @brief Parse TOP file
     * 
     * @param filename Path to TOP file
     * @return Topology Parsed topology information
     */
    static Topology parse(const std::string& filename);

private:
    /**
     * @brief Parse [ atomtypes ] section
     * 
     * @param is Input stream positioned at the start of atomtypes section
     * @param top Topology object to store parsed information
     * @return true if parsing successful
     * @return false if parsing failed
     */
    static bool parse_atomtypes_section(std::istream& is, Topology& top);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_TOP_PARSER_HPP
