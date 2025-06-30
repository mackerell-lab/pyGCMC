// src/io/topology/topParserMoleculeBuilder.cpp

#include "topParserMoleculeBuilder.hpp"
#include "topParserSections.hpp"
#include "topParserUtilities.hpp"

namespace pygcmc {
namespace io {

bool TopParserMoleculeBuilder::build_molecules(model::Topology& topology, const SectionProcessorResult& result) {
    TopParserUtilities::debug_print("\n=== Counting atoms per molecule type ===\n");
    std::map<std::string, int> atoms_per_molecule = count_atoms_per_molecule(result);

    TopParserUtilities::debug_print("\n=== Adding molecules to topology ===\n");
    int atom_offset = 0;
    int mol_index = 0;

    // Add molecules in order according to [ molecules ] section
    for (const auto& mol : result.molecule_order) {
        const std::string& mol_type = mol.first;
        int count = mol.second;

        // Skip if we don't have the molecule definition
        if (result.molecule_atoms_temp.find(mol_type) == result.molecule_atoms_temp.end()) {
            TopParserUtilities::debug_print("Warning: No definition found for molecule type ", mol_type, "\n");
            continue;
        }

        if (!add_molecule_copies(topology, mol_type, count, mol_index, atom_offset, 
                               atoms_per_molecule, result)) {
            return false;
        }
    }

    TopParserUtilities::debug_print("\n=== Final topology statistics ===\n");
    TopParserUtilities::debug_print("Total atoms: ", topology.get_num_atoms(), "\n");
    TopParserUtilities::debug_print("Total bonds: ", topology.get_num_bonds(), "\n");
    TopParserUtilities::debug_print("Total angles: ", topology.get_num_angles(), "\n");
    TopParserUtilities::debug_print("Total dihedrals: ", topology.get_num_dihedrals(), "\n");
    TopParserUtilities::debug_print("Total impropers: ", topology.get_num_impropers(), "\n");
    TopParserUtilities::debug_print("Total residues: ", topology.get_num_residues(), "\n");
    TopParserUtilities::debug_print("Total segments: ", topology.get_num_segments(), "\n");

    return true;
}

std::map<std::string, int> TopParserMoleculeBuilder::count_atoms_per_molecule(const SectionProcessorResult& result) {
    std::map<std::string, int> atoms_per_molecule;

    // Count atoms per molecule type, but only for molecules we'll actually use
    for (const auto& mol_type : result.used_molecule_types) {
        auto it = result.molecule_atoms_temp.find(mol_type);
        if (it != result.molecule_atoms_temp.end()) {
            int atom_count = 0;
            TopParserUtilities::debug_print("Checking atoms for molecule type ", mol_type, ":\n");
            for (const auto& line_info : it->second) {
                auto tokens = TopParserUtilities::split(TopParserUtilities::remove_comment(line_info.content));
                if (tokens.size() >= 8) {
                    atom_count++;
                    TopParserUtilities::debug_print("  Atom ", tokens[0], " (", tokens[1], 
                             ") in residue ", tokens[3], "\n");
                }
            }
            atoms_per_molecule[mol_type] = atom_count;
            TopParserUtilities::debug_print("Molecule ", mol_type, " has ", atom_count, " atoms\n");
        } else {
            TopParserUtilities::debug_print("Warning: No atom definitions found for molecule type ", mol_type, "\n");
        }
    }

    return atoms_per_molecule;
}

bool TopParserMoleculeBuilder::add_molecule_copies(model::Topology& topology,
                                                 const std::string& mol_type,
                                                 int count,
                                                 int& mol_index,
                                                 int& atom_offset,
                                                 const std::map<std::string, int>& atoms_per_molecule,
                                                 const SectionProcessorResult& result) {
    TopParserUtilities::debug_print("Adding ", count, " copies of molecule ", mol_type, "\n");

    // Add 'count' copies of this molecule type
    for (int i = 0; i < count; i++) {
        // Create unique segment name for this molecule instance
        std::string segment_name = mol_type + "_" + std::to_string(mol_index++);
        topology.add_segment(segment_name);
        TopParserUtilities::debug_print("  Creating segment ", segment_name, "\n");

        if (!add_molecule_instance(topology, mol_type, segment_name, atom_offset, result)) {
            return false;
        }

        // Update atom offset for next molecule instance
        auto it = atoms_per_molecule.find(mol_type);
        if (it != atoms_per_molecule.end()) {
            atom_offset += it->second;
        }
    }

    return true;
}

bool TopParserMoleculeBuilder::add_molecule_instance(model::Topology& topology,
                                                   const std::string& mol_type,
                                                   const std::string& segment_name,
                                                   int atom_offset,
                                                   const SectionProcessorResult& result) {
    // Add atoms for this molecule instance
    auto atoms_it = result.molecule_atoms_temp.find(mol_type);
    if (atoms_it != result.molecule_atoms_temp.end()) {
        if (!TopParserSections::parse_atoms_section(atoms_it->second, topology, segment_name)) {
            const auto& line = atoms_it->second.front();
            TopParserUtilities::debug_print("Error parsing atoms for molecule ", mol_type, 
                     " at ", line.source_file, ":", line.line_number, "\n");
            return false;
        }
        TopParserUtilities::debug_print("  Added atoms for segment ", segment_name, "\n");
    }

    // Add bonds for this molecule instance
    auto bonds_it = result.molecule_bonds_temp.find(mol_type);
    if (bonds_it != result.molecule_bonds_temp.end()) {
        if (!TopParserSections::parse_bonds_section(bonds_it->second, topology, atom_offset)) {
            const auto& line = bonds_it->second.front();
            TopParserUtilities::debug_print("Error parsing bonds for molecule ", mol_type, 
                     " at ", line.source_file, ":", line.line_number, "\n");
            return false;
        }
        TopParserUtilities::debug_print("  Added bonds for segment ", segment_name, "\n");
    }

    // Add angles for this molecule instance
    auto angles_it = result.molecule_angles_temp.find(mol_type);
    if (angles_it != result.molecule_angles_temp.end()) {
        if (!TopParserSections::parse_angles_section(angles_it->second, topology, atom_offset)) {
            const auto& line = angles_it->second.front();
            TopParserUtilities::debug_print("Error parsing angles for molecule ", mol_type, 
                     " at ", line.source_file, ":", line.line_number, "\n");
            return false;
        }
        TopParserUtilities::debug_print("  Added angles for segment ", segment_name, "\n");
    }

    // Add dihedrals for this molecule instance
    auto dihedrals_it = result.molecule_dihedrals_temp.find(mol_type);
    if (dihedrals_it != result.molecule_dihedrals_temp.end()) {
        if (!TopParserSections::parse_dihedrals_section(dihedrals_it->second, topology, atom_offset)) {
            const auto& line = dihedrals_it->second.front();
            TopParserUtilities::debug_print("Error parsing dihedrals for molecule ", mol_type, 
                     " at ", line.source_file, ":", line.line_number, "\n");
            return false;
        }
        TopParserUtilities::debug_print("  Added dihedrals for segment ", segment_name, "\n");
    }

    // Add impropers for this molecule instance
    auto impropers_it = result.molecule_impropers_temp.find(mol_type);
    if (impropers_it != result.molecule_impropers_temp.end()) {
        if (!TopParserSections::parse_impropers_section(impropers_it->second, topology, atom_offset)) {
            const auto& line = impropers_it->second.front();
            TopParserUtilities::debug_print("Error parsing impropers for molecule ", mol_type, 
                     " at ", line.source_file, ":", line.line_number, "\n");
            return false;
        }
        TopParserUtilities::debug_print("  Added impropers for segment ", segment_name, "\n");
    }

    // Add CMAPs for this molecule instance
    auto cmaps_it = result.molecule_cmaps_temp.find(mol_type);
    if (cmaps_it != result.molecule_cmaps_temp.end()) {
        if (!TopParserSections::parse_cmaps_section(cmaps_it->second, topology, atom_offset)) {
            const auto& line = cmaps_it->second.front();
            TopParserUtilities::debug_print("Error parsing CMAPs for molecule ", mol_type, 
                     " at ", line.source_file, ":", line.line_number, "\n");
            return false;
        }
        TopParserUtilities::debug_print("  Added CMAPs for segment ", segment_name, "\n");
    }

    return true;
}

} // namespace io
} // namespace pygcmc