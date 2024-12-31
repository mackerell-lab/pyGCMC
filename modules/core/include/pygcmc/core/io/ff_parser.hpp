// modules/core/include/pygcmc/core/io/ff_parser.hpp

#ifndef PYGCMC_CORE_IO_FF_PARSER_HPP
#define PYGCMC_CORE_IO_FF_PARSER_HPP

#include <string>
#include <unordered_map>
#include <utility>
#include <stdexcept>

namespace pygcmc {
namespace core {
namespace io {

// Structure to store Lennard-Jones parameters
struct ForceFieldPair {
    double epsilon;  // well depth
    double rmin;     // Rmin/2 in CHARMM format

    ForceFieldPair(double e=0.0, double r=0.0) : epsilon(e), rmin(r) {}
};

// Hash function for pair of strings (atom types)
struct PairStringHash {
    std::size_t operator()(const std::pair<std::string, std::string>& p) const {
        auto h1 = std::hash<std::string>()(p.first);
        auto h2 = std::hash<std::string>()(p.second);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

class FFParser {
public:
    // Parse CHARMM parameter file
    bool parse(const std::string& filename);

    // Get nonbonded parameters
    const std::unordered_map<std::string, ForceFieldPair>& get_nonbonded_params() const {
        return nonbonded_params_;
    }

    // Get NBFIX parameters
    const std::unordered_map<
        std::pair<std::string, std::string>, 
        ForceFieldPair, 
        PairStringHash
    >& get_nbfix_params() const {
        return nbfix_params_;
    }

private:
    // Storage for parameters
    std::unordered_map<std::string, ForceFieldPair> nonbonded_params_;
    std::unordered_map<
        std::pair<std::string, std::string>, 
        ForceFieldPair, 
        PairStringHash
    > nbfix_params_;

    // Helper functions for parsing sections
    void parse_nonbonded_section(std::istream& in);
    void parse_nbfix_section(std::istream& in);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_FF_PARSER_HPP

