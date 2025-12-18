#include "ItpNonbondedParser.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace pygcmc {
namespace io {

namespace {

std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        ++start;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }
    return s.substr(start, end - start);
}

std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool startsWith(const std::string& s, const char* prefix) {
    const size_t n = std::char_traits<char>::length(prefix);
    return s.size() >= n && s.compare(0, n, prefix) == 0;
}

std::string stripComment(const std::string& line) {
    const size_t pos = line.find(';');
    if (pos == std::string::npos) {
        return line;
    }
    return line.substr(0, pos);
}

std::vector<std::string> splitWhitespace(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string t;
    while (iss >> t) {
        tokens.push_back(t);
    }
    return tokens;
}

std::pair<std::string, std::string> canonicalPair(std::string a, std::string b) {
    if (b < a) {
        std::swap(a, b);
    }
    return {a, b};
}

bool isSectionHeader(const std::string& line) {
    return !line.empty() && line.front() == '[' && line.find(']') != std::string::npos;
}

std::string parseSectionName(const std::string& headerLine) {
    const auto close = headerLine.find(']');
    std::string inside = headerLine.substr(1, close - 1);
    return toLower(trim(inside));
}

} // namespace

ItpNonbondedParser::Result ItpNonbondedParser::parse_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open ITP file: " + filename);
    }

    Result result;
    // Track overrides by section so we can apply a deterministic precedence.
    // For compatibility with the planned "strict-gromacs" path (and to match the intended
    // gcmc_gpu behavior), treat nonbond_params as NBFIX and let it override pairtypes.
    std::map<std::pair<std::string, std::string>, LJ> pairtypesOverrides;
    std::map<std::pair<std::string, std::string>, LJ> nbfixOverrides;
    std::string section;
    std::string line;

    while (std::getline(file, line)) {
        line = trim(stripComment(line));
        if (line.empty()) {
            continue;
        }

        if (isSectionHeader(line)) {
            section = parseSectionName(line);
            continue;
        }

        // gcmc_gpu-compatible preprocessor handling:
        // - When encountering "#ifdef", skip the first branch and parse the #else branch if present.
        if (startsWith(line, "#ifdef")) {
            std::string directiveLine;
            while (std::getline(file, directiveLine)) {
                directiveLine = trim(stripComment(directiveLine));
                if (startsWith(directiveLine, "#else") || startsWith(directiveLine, "#endif")) {
                    break;
                }
            }
            // Skip directive line itself; if it was "#else", parsing resumes in the else branch.
            continue;
        }
        if (startsWith(line, "#else") || startsWith(line, "#endif") || startsWith(line, "#include")) {
            continue;
        }

        const auto tokens = splitWhitespace(line);
        if (tokens.empty()) {
            continue;
        }

        if (section == "atomtypes") {
            // Typical format: name bond_type mass charge ptype sigma epsilon
            if (tokens.size() < 3) {
                continue;
            }
            const std::string& typeName = tokens[0];
            try {
                const double sigma = std::stod(tokens[tokens.size() - 2]);
                const double eps = std::stod(tokens[tokens.size() - 1]);
                result.atomTypes[typeName] = LJ{sigma, eps};
            } catch (const std::exception&) {
                continue;
            }
        } else if (section == "pairtypes" || section == "nonbond_params") {
            // Typical format: type1 type2 func sigma epsilon
            if (tokens.size() < 4) {
                continue;
            }
            const std::string& t1 = tokens[0];
            const std::string& t2 = tokens[1];
            try {
                const double sigma = std::stod(tokens[tokens.size() - 2]);
                const double eps = std::stod(tokens[tokens.size() - 1]);
                const auto key = canonicalPair(t1, t2);
                if (section == "pairtypes") {
                    pairtypesOverrides[key] = LJ{sigma, eps};
                } else {
                    nbfixOverrides[key] = LJ{sigma, eps};
                }
            } catch (const std::exception&) {
                continue;
            }
        }
    }

    // Merge overrides with precedence: pairtypes first, then nonbond_params (NBFIX).
    result.pairOverrides = std::move(pairtypesOverrides);
    for (const auto& [pair, lj] : nbfixOverrides) {
        result.pairOverrides[pair] = lj;
    }

    return result;
}

ItpNonbondedParser::Result ItpNonbondedParser::parse_files(const std::vector<std::string>& filenames) {
    Result merged;
    for (const auto& f : filenames) {
        Result r = parse_file(f);
        for (const auto& [name, lj] : r.atomTypes) {
            merged.atomTypes[name] = lj;
        }
        for (const auto& [pair, lj] : r.pairOverrides) {
            merged.pairOverrides[pair] = lj;
        }
    }
    return merged;
}

} // namespace io
} // namespace pygcmc
