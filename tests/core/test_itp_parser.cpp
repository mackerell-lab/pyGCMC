// tests/core/test_itp_parser.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/io/itp_parser.hpp"
#include "config.hpp"
#include <cmath>

using namespace pygcmc::core::io;

TEST(ITPParserTest, ParseBENXFile) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/benx.itp";
    ASSERT_TRUE(parser.parse(itp_file));

    // Test getting atom properties for BENX residue
    double charge, mass;
    
    // Test CG atom
    ASSERT_TRUE(parser.get_atom_properties("BENX", "CG", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.115);
    EXPECT_DOUBLE_EQ(mass, 12.011);

    // Test HG atom
    ASSERT_TRUE(parser.get_atom_properties("BENX", "HG", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.115);
    EXPECT_DOUBLE_EQ(mass, 1.008);

    // Test CD1 atom
    ASSERT_TRUE(parser.get_atom_properties("BENX", "CD1", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.115);
    EXPECT_DOUBLE_EQ(mass, 12.011);

    // Test non-existent atom
    ASSERT_FALSE(parser.get_atom_properties("BENX", "XXX", charge, mass));
}

TEST(ITPParserTest, UpdateBENXAtoms) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/benx.itp";
    ASSERT_TRUE(parser.parse(itp_file));

    // Create test PDB atoms for BENX
    std::vector<PDBAtom> pdb_atoms;
    
    // Add BENX atoms
    PDBAtom cg(1, "CG", "BENX", 1);
    PDBAtom hg(2, "HG", "BENX", 1);
    PDBAtom cd1(3, "CD1", "BENX", 1);
    
    pdb_atoms.push_back(cg);
    pdb_atoms.push_back(hg);
    pdb_atoms.push_back(cd1);

    // Update atoms
    int updated = parser.update_pdb_atoms(pdb_atoms);
    EXPECT_EQ(updated, 3);

    // Verify updated atoms
    EXPECT_EQ(pdb_atoms[0].topo_type, "CG2R61");
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_charge, -0.115);
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_mass, 12.011);

    EXPECT_EQ(pdb_atoms[1].topo_type, "HGR61");
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_charge, 0.115);
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_mass, 1.008);

    EXPECT_EQ(pdb_atoms[2].topo_type, "CG2R61");
    EXPECT_DOUBLE_EQ(pdb_atoms[2].topo_charge, -0.115);
    EXPECT_DOUBLE_EQ(pdb_atoms[2].topo_mass, 12.011);
}

TEST(ITPParserTest, ParseSOLFile) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/sol.itp";
    ASSERT_TRUE(parser.parse(itp_file));

    // Test getting atom properties for SOL residue
    double charge, mass;
    
    // Test OW atom
    ASSERT_TRUE(parser.get_atom_properties("SOL", "OW", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.834);
    EXPECT_DOUBLE_EQ(mass, 15.9994);

    // Test HW1 atom
    ASSERT_TRUE(parser.get_atom_properties("SOL", "HW1", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.417);
    EXPECT_DOUBLE_EQ(mass, 1.008);

    // Test HW2 atom
    ASSERT_TRUE(parser.get_atom_properties("SOL", "HW2", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.417);
    EXPECT_DOUBLE_EQ(mass, 1.008);
}

TEST(ITPParserTest, ParsePRPFile) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/prpx.itp";
    ASSERT_TRUE(parser.parse(itp_file));

    // Test getting atom properties for PRP residue
    double charge, mass;
    
    // Test C1 atom
    ASSERT_TRUE(parser.get_atom_properties("PRPX", "C1", charge, mass));
    EXPECT_DOUBLE_EQ(charge, -0.27);
    EXPECT_DOUBLE_EQ(mass, 12.011);

    // Test H11 atom
    ASSERT_TRUE(parser.get_atom_properties("PRPX", "H11", charge, mass));
    EXPECT_DOUBLE_EQ(charge, 0.09);
    EXPECT_DOUBLE_EQ(mass, 1.008);
}

TEST(ITPParserTest, MissingTopologyInfo) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/benx.itp";
    ASSERT_TRUE(parser.parse(itp_file));

    std::vector<PDBAtom> pdb_atoms;
    
    // Add existing atom
    PDBAtom cg(1, "CG", "BENX", 1);
    // Add non-existent atom
    PDBAtom xxx(2, "XXX", "BENX", 1);
    // Add atom from different residue
    PDBAtom ow(3, "OW", "SOL", 1);
    
    pdb_atoms.push_back(cg);
    pdb_atoms.push_back(xxx);
    pdb_atoms.push_back(ow);

    auto missing_info = parser.get_missing_topology_info(pdb_atoms);
    
    // BENX XXX should be missing
    EXPECT_TRUE(missing_info["BENX"].find("XXX") != missing_info["BENX"].end());
    // SOL OW should be missing (different residue)
    EXPECT_TRUE(missing_info["SOL"].find("OW") != missing_info["SOL"].end());
    // BENX CG should not be missing
    EXPECT_TRUE(missing_info["BENX"].find("CG") == missing_info["BENX"].end());
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}