// src/io/topology/psfParserSections.cpp

#include "psfParserSections.hpp"
#include "psfParserStringUtils.hpp"
#include <sstream>
#include <iostream>
#include <array>
#include <algorithm>

namespace {
    // Local trim function to match original behavior exactly
    std::string trim(const std::string& str) {
        const auto start = str.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) {
            return "";
        }
        const auto end = str.find_last_not_of(" \t\r\n");
        return str.substr(start, end - start + 1);
    }
}

namespace pygcmc {
namespace io {

bool PSFParserSections::parse_impropers_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    int num_impropers = 0;
    if (!(iss >> num_impropers)) {
        std::cerr << "Failed to parse number of impropers" << std::endl;
        return false;
    }

    if (num_impropers <= 0) {
        // no impropers
        return true;
    }

    int impropers_added = 0;
    std::vector<int> tmp_indices;

    for (size_t i = 1; i < lines.size(); ++i) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) continue;
            idx--; // 0-based
            if (idx < 0 || idx >= topology.get_num_atoms()) {
                std::cerr << "Invalid atom index in impropers" << std::endl;
                return false;
            }
            tmp_indices.push_back(idx);
            if (tmp_indices.size() == 4) {
                topology.add_improper(tmp_indices[0],
                                    tmp_indices[1],
                                    tmp_indices[2],
                                    tmp_indices[3]);
                tmp_indices.clear();
                impropers_added++;
            }
        }
    }

    return impropers_added == num_impropers;
}

bool PSFParserSections::parse_donors_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_donors = 0;
    if (!(iss >> num_donors)) {
        std::cerr << "Failed to parse number of donors" << std::endl;
        return false;
    }

    size_t donors_parsed = 0;
    for (size_t i = 1; i < lines.size() && donors_parsed < num_donors; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int donor_idx, hydrogen_idx;
        while (iss_line >> donor_idx >> hydrogen_idx) {
            if (donor_idx == 0 || hydrogen_idx == 0) continue;
            donor_idx--; // Convert to 0-based indexing
            hydrogen_idx--;
            topology.add_donor(donor_idx, hydrogen_idx);
            donors_parsed++;
            if (donors_parsed == num_donors) {
                break;
            }
        }
    }
    return donors_parsed == num_donors;
}

bool PSFParserSections::parse_acceptors_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_acceptors = 0;
    if (!(iss >> num_acceptors)) {
        std::cerr << "Failed to parse number of acceptors" << std::endl;
        return false;
    }

    size_t acceptors_parsed = 0;
    for (size_t i = 1; i < lines.size() && acceptors_parsed < num_acceptors; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int acceptor_idx;
        while (iss_line >> acceptor_idx) {
            if (acceptor_idx == 0) continue;
            acceptor_idx--; // Convert to 0-based indexing
            topology.add_acceptor(acceptor_idx);
            acceptors_parsed++;
            if (acceptors_parsed == num_acceptors) {
                break;
            }
        }
    }
    return acceptors_parsed == num_acceptors;
}

bool PSFParserSections::parse_cmap_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_cmaps = 0;
    if (!(iss >> num_cmaps)) {
        std::cerr << "Failed to parse number of CMAP terms" << std::endl;
        return false;
    }

    size_t cmaps_parsed = 0;
    std::vector<int> buffer;
    for (size_t i = 1; i < lines.size() && cmaps_parsed < num_cmaps; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) continue; // Skip placeholder zeros
            buffer.push_back(idx - 1); // Convert to 0-based indexing
            if (buffer.size() == 8) {
                std::array<int, 8> cmap_array;
                std::copy_n(buffer.begin(), 8, cmap_array.begin());
                topology.add_cmap(cmap_array);
                buffer.clear();
                cmaps_parsed++;
                if (cmaps_parsed == num_cmaps) {
                    break;
                }
            }
        }
    }
    return cmaps_parsed == num_cmaps;
}

bool PSFParserSections::parse_groups_from_lines(const std::vector<std::string>& lines, model::Topology& topology) {
    std::istringstream iss(trim(lines[0]));
    size_t num_groups = 0;
    if (!(iss >> num_groups)) {
        std::cerr << "Failed to parse number of groups" << std::endl;
        return false;
    }

    size_t groups_parsed = 0;
    std::vector<int> current_group;
    int current_group_id = 1;  // Start with group ID 1

    for (size_t i = 1; i < lines.size() && groups_parsed < num_groups; i++) {
        std::string line = trim(lines[i]);
        std::istringstream iss_line(line);
        int idx;
        while (iss_line >> idx) {
            if (idx == 0) {
                // Zero marks end of current group
                if (!current_group.empty()) {
                    topology.add_group(current_group_id++, current_group, "");  // Empty type string
                    current_group.clear();
                    groups_parsed++;
                    if (groups_parsed == num_groups) {
                        break;
                    }
                }
                continue;
            }
            current_group.push_back(idx - 1); // Convert to 0-based indexing
        }
    }
    // Add last group if not empty
    if (!current_group.empty()) {
        topology.add_group(current_group_id, current_group, "");
        groups_parsed++;
    }
    return groups_parsed == num_groups;
}

} // namespace io
} // namespace pygcmc