#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace pygcmc {
namespace io {

/**
 * @brief Minimal GROMACS .itp nonbonded parser used by gcmc_gpu-style INP.
 *
 * Supported sections:
 * - [ atomtypes ]       -> per-type sigma/epsilon (nm, kJ/mol)
 * - [ pairtypes ]       -> pair overrides (sigma/epsilon)
 * - [ nonbond_params ]  -> pair overrides (sigma/epsilon)
 *
 * Notes:
 * - This parser intentionally focuses on the columns used by tmp/gcmc_gpu/source/parse.cpp.
 * - It does not currently evaluate preprocessor symbols; #ifdef blocks are handled in a
 *   gcmc_gpu-compatible way (skip the first branch, parse the #else branch if present).
 */
class ItpNonbondedParser {
public:
    struct LJ {
        double sigma_nm{0.0};
        double epsilon_kj{0.0};
    };

    struct Result {
        std::map<std::string, LJ> atomTypes;
        std::map<std::pair<std::string, std::string>, LJ> pairOverrides;
    };

    static Result parse_file(const std::string& filename);
    static Result parse_files(const std::vector<std::string>& filenames);
};

} // namespace io
} // namespace pygcmc

