// src/io/topology/TopParserStringUtils.hpp

#pragma once

#include "../../model/ModelModule.hpp"
#include <string>

namespace pygcmc {
namespace io {

class TopParserStringUtils {
public:
    /**
     * Parse a topology string using an isolated temporary file
     * @param top_str The topology string content
     * @return Parsed topology
     */
    static model::Topology parse_string(const std::string& top_str);
};

} // namespace io
} // namespace pygcmc
