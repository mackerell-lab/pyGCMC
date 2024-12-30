// tests/core/test_top_parser.cpp

#include <gtest/gtest.h>
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "config.hpp"

using namespace pygcmc::core::io;

TEST(TopParserTest, ParseTopologyFile) {
    TopParser parser;
    std::string top_file = std::string(PDB_DATA_DIR) + "/test.top";
    ASSERT_TRUE(parser.parse(top_file));

    // Test getting atom properties for ALA residue
    double charge, mass;
    // Note: The last defined ALA N atom (residue 10) will be used
    ASSERT_TRUE(parser.get_atom_properties("ALA", "N", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.47);  // From residue 10
    EXPECT_DOUBLE_EQ(mass, 14.007);

    // Test VAL residue atoms
    ASSERT_TRUE(parser.get_atom_properties("VAL", "CA", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.07);
    EXPECT_DOUBLE_EQ(mass, 12.011);

    ASSERT_TRUE(parser.get_atom_properties("VAL", "CG1", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.27);
    EXPECT_DOUBLE_EQ(mass, 12.011);

    // Test PRO residue atoms
    ASSERT_TRUE(parser.get_atom_properties("PRO", "N", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.29);
    EXPECT_DOUBLE_EQ(mass, 14.007);

    ASSERT_TRUE(parser.get_atom_properties("PRO", "CD", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.0);
    EXPECT_DOUBLE_EQ(mass, 12.011);

    // Test ASN residue atoms
    ASSERT_TRUE(parser.get_atom_properties("ASN", "CG", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.55);
    EXPECT_DOUBLE_EQ(mass, 12.011);

    ASSERT_TRUE(parser.get_atom_properties("ASN", "OD1", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.55);
    EXPECT_DOUBLE_EQ(mass, 15.9994);

    // Test GLN residue atoms
    ASSERT_TRUE(parser.get_atom_properties("GLN", "NE2", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.62);
    EXPECT_DOUBLE_EQ(mass, 14.007);

    ASSERT_TRUE(parser.get_atom_properties("GLN", "HE21", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.32);
    EXPECT_DOUBLE_EQ(mass, 1.008);

    // Test non-existent atoms
    ASSERT_FALSE(parser.get_atom_properties("XXX", "XXX", charge, mass));
    ASSERT_FALSE(parser.get_atom_properties("BEN", "C1", charge, mass));
    ASSERT_FALSE(parser.get_atom_properties("SOL", "OW", charge, mass));
}

TEST(TopParserTest, UpdatePDBAtoms) {
    // First parse the PDB file
    std::string pdb_file = std::string(PDB_DATA_DIR) + "/test.pdb";
    auto [cell_params, residues] = PDBParser::parse(pdb_file);
    ASSERT_FALSE(residues.empty());

    // Extract atoms from residues
    std::vector<PDBAtom> pdb_atoms;
    for (const auto& residue : residues) {
        pdb_atoms.insert(pdb_atoms.end(), residue.atoms.begin(), residue.atoms.end());
    }

    // Then parse the topology file and update PDB atoms
    TopParser top_parser;
    std::string top_file = std::string(PDB_DATA_DIR) + "/test.top";
    ASSERT_TRUE(top_parser.parse(top_file));
    
    int updated = top_parser.update_pdb_atoms(pdb_atoms);
    EXPECT_GT(updated, 0);

    // Verify updated atoms
    for (const auto& atom : pdb_atoms) {
        if (atom.residue == "ALA" && atom.name == "N") {
            if (atom.sequence == 7) {  // First ALA residue
                EXPECT_DOUBLE_EQ(atom.topo_charge, -0.3);
                EXPECT_DOUBLE_EQ(atom.topo_mass, 14.007);
                EXPECT_EQ(atom.topo_type, "NH3");
            } else if (atom.sequence == 10) {  // Second ALA residue
                EXPECT_DOUBLE_EQ(atom.topo_charge, -0.47);
                EXPECT_DOUBLE_EQ(atom.topo_mass, 14.007);
                EXPECT_EQ(atom.topo_type, "NH1");
            }
        }
        else if (atom.residue == "VAL" && atom.name == "CA") {
            EXPECT_DOUBLE_EQ(atom.topo_charge, 0.07);
            EXPECT_DOUBLE_EQ(atom.topo_mass, 12.011);
            EXPECT_EQ(atom.topo_type, "CT1");
        }
        else if (atom.residue == "PRO" && atom.name == "N") {
            EXPECT_DOUBLE_EQ(atom.topo_charge, -0.29);
            EXPECT_DOUBLE_EQ(atom.topo_mass, 14.007);
            EXPECT_EQ(atom.topo_type, "N");
        }
        else if (atom.residue == "ASN" && atom.name == "CG") {
            EXPECT_DOUBLE_EQ(atom.topo_charge, 0.55);
            EXPECT_DOUBLE_EQ(atom.topo_mass, 12.011);
            EXPECT_EQ(atom.topo_type, "CC");
        }
        else if (atom.residue == "GLN" && atom.name == "NE2") {
            EXPECT_DOUBLE_EQ(atom.topo_charge, -0.62);
            EXPECT_DOUBLE_EQ(atom.topo_mass, 14.007);
            EXPECT_EQ(atom.topo_type, "NH2");
        }
    }
} 