// src/io/topology/topParserSectionProcessor.hpp

#pragma once

#include "topParserStructures.hpp"
#include <string>
#include <vector>
#include <map>
#include <set>

namespace pygcmc {
namespace io {

struct SectionProcessorResult {
    std::vector<std::string> molecule_types_order;
    std::map<std::string, std::vector<LineInfo>> molecule_atoms_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_bonds_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_angles_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_dihedrals_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_impropers_temp;
    std::map<std::string, std::vector<LineInfo>> molecule_cmaps_temp;
    std::vector<LineInfo> molecules_lines;
    std::vector<std::pair<std::string, int>> molecule_order;
    std::set<std::string> used_molecule_types;
    bool found_valid_section = false;
};

class TopParserSectionProcessor {
public:
    /**
     * Process all sections from collected lines
     * @param all_lines All preprocessed lines from the topology file
     * @return Processed section data
     */
    static SectionProcessorResult process_sections(const std::vector<LineInfo>& all_lines);

private:
    static void process_moleculetype_section(const std::string& trimmed, 
                                           [[maybe_unused]] const LineInfo& line_info,
                                           std::string& current_mol_type,
                                           std::vector<std::string>& molecule_types_order,
                                           std::map<std::string, std::vector<LineInfo>>& molecule_atoms_temp,
                                           std::map<std::string, std::vector<LineInfo>>& molecule_bonds_temp,
                                           std::map<std::string, std::vector<LineInfo>>& molecule_angles_temp,
                                           std::map<std::string, std::vector<LineInfo>>& molecule_dihedrals_temp,
                                           std::map<std::string, std::vector<LineInfo>>& molecule_impropers_temp);

    static void process_molecules_section(const std::vector<LineInfo>& molecules_lines,
                                        const std::vector<std::string>& molecule_types_order,
                                        std::vector<std::pair<std::string, int>>& molecule_order,
                                        std::set<std::string>& used_molecule_types);
};

} // namespace io
} // namespace pygcmc