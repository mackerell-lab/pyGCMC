// src/io/topology/topParserMoleculeBuilder.hpp

#pragma once

#include "../../model/ModelModule.hpp"
#include "topParserSectionProcessor.hpp"
#include <string>
#include <vector>
#include <map>
#include <set>

namespace pygcmc {
namespace io {

class TopParserMoleculeBuilder {
public:
    /**
     * Build molecules and add them to topology based on processed section data
     * @param topology The topology to add molecules to
     * @param result The processed section data
     * @return True if successful, false otherwise
     */
    static bool build_molecules(model::Topology& topology, const SectionProcessorResult& result);

private:
    static std::map<std::string, int> count_atoms_per_molecule(const SectionProcessorResult& result);
    
    static bool add_molecule_copies(model::Topology& topology,
                                  const std::string& mol_type,
                                  int count,
                                  int& mol_index,
                                  int& atom_offset,
                                  const std::map<std::string, int>& atoms_per_molecule,
                                  const SectionProcessorResult& result);
    
    static bool add_molecule_instance(model::Topology& topology,
                                    const std::string& mol_type,
                                    const std::string& segment_name,
                                    int atom_offset,
                                    const SectionProcessorResult& result);
};

} // namespace io
} // namespace pygcmc