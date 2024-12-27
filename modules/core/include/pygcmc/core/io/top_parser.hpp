// modules/core/include/pygcmc/core/io/top_parser.hpp

#ifndef PYGCMC_CORE_IO_TOP_PARSER_HPP
#define PYGCMC_CORE_IO_TOP_PARSER_HPP

#include "parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief TOP Parser
 * 
 * Parses TOP files to extract atom type information and force field parameters.
 */
class TopParser {
public:
    /**
     * @brief Parse a TOP file
     * 
     * @param filename Path to the TOP file
     * @return Topology Parsed topology information
     */
    static Topology parse(const std::string& filename);

private:
    /**
     * @brief Parse [ atomtypes ] section
     * 
     * @param line Line from the TOP file
     * @param top Topology object to store parsed information
     * @return true if parsing successful
     * @return false if parsing failed
     */
    static bool parse_atomtypes_section(const std::string& line, Topology& top);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_TOP_PARSER_HPP
