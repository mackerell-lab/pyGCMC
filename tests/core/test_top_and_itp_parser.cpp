// tests/core/test_top_and_itp_parser.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/io/top_parser.hpp"
#include "config.hpp"
#include <filesystem>

using namespace pygcmc::core::io;

class TopParserTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Get the test data directory from environment or use default
        const char* test_data_dir = std::getenv("TEST_DATA_DIR");
        if (test_data_dir) {
            data_dir_ = test_data_dir;
        } else {
            data_dir_ = std::string(PDB_DATA_DIR);
        }
    }

    std::string data_dir_;
};

TEST_F(TopParserTest, ParseTopWithoutIncludes) {
    TopParser parser;
    std::string top_file = data_dir_ + "/test.top";
    ASSERT_TRUE(parser.parse(top_file)) << "Failed to parse topology file";
    
    // Verify some basic atom properties from the main topology
    double charge, mass;
    EXPECT_TRUE(parser.get_atom_properties("ALA", "N", charge, mass));
    EXPECT_NEAR(charge, -0.3, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);

    // SOL should not be found when includes are disabled
    EXPECT_FALSE(parser.get_atom_properties("SOL", "OW", charge, mass));
}

TEST_F(TopParserTest, ParseTopWithIncludes) {
    TopParser parser;
    std::string top_file = data_dir_ + "/test.top";
    ASSERT_TRUE(parser.parse_with_includes(top_file)) << "Failed to parse topology file with includes";
    
    // Verify atoms from main topology
    double charge, mass;
    EXPECT_TRUE(parser.get_atom_properties("ALA", "N", charge, mass));
    EXPECT_NEAR(charge, -0.3, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);
    
    // Verify atoms from included files
    EXPECT_TRUE(parser.get_atom_properties("SOL", "OW", charge, mass));
    EXPECT_NEAR(charge, -0.834, 1e-6);
    EXPECT_NEAR(mass, 15.9994, 1e-6);
}

TEST_F(TopParserTest, HandleMissingIncludes) {
    TopParser parser;
    std::string top_file = data_dir_ + "/test_missing_includes.top";
    
    // Should not fail even with missing includes
    EXPECT_TRUE(parser.parse_with_includes(top_file)) << "Parser failed with missing includes";
    
    // Should still be able to read atoms from main file
    double charge, mass;
    EXPECT_TRUE(parser.get_atom_properties("ALA", "N", charge, mass));
    EXPECT_NEAR(charge, -0.3, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);
}

TEST_F(TopParserTest, UpdatePDBAtoms) {
    TopParser parser;
    std::string top_file = data_dir_ + "/test.top";
    ASSERT_TRUE(parser.parse_with_includes(top_file)) << "Failed to parse topology file";
    
    std::vector<PDBAtom> atoms;
    PDBAtom atom;
    atom.residue = "SOL";
    atom.name = "OW";
    atom.sequence = 1;
    atoms.push_back(atom);
    
    EXPECT_EQ(parser.update_pdb_atoms(atoms), 1) << "Failed to update PDB atoms";
    EXPECT_NEAR(atoms[0].topo_charge, -0.834, 1e-6);
    EXPECT_NEAR(atoms[0].topo_mass, 15.9994, 1e-6);
}

TEST_F(TopParserTest, GetMissingTopologyInfo) {
    TopParser parser;
    std::string top_file = data_dir_ + "/test.top";
    ASSERT_TRUE(parser.parse_with_includes(top_file)) << "Failed to parse topology file";
    
    std::vector<PDBAtom> atoms;
    PDBAtom atom;
    atom.residue = "UNKNOWN";
    atom.name = "X";
    atom.sequence = 1;
    atoms.push_back(atom);
    
    auto missing = parser.get_missing_topology_info(atoms);
    EXPECT_EQ(missing.size(), 1) << "Expected one missing residue";
    EXPECT_EQ(missing["UNKNOWN"].size(), 1) << "Expected one missing atom";
    EXPECT_TRUE(missing["UNKNOWN"].find("X") != missing["UNKNOWN"].end());
}

TEST_F(TopParserTest, ParseTopologyWithIncludes) {
    TopParser parser;
    std::string top_file = data_dir_ + "/test.top";
    ASSERT_TRUE(parser.parse_with_includes(top_file)) << "Failed to parse topology file with includes";
    
    double charge, mass;
    // Test atoms from main file
    ASSERT_TRUE(parser.get_atom_properties("ALA", "N", charge, mass));
    EXPECT_NEAR(charge, -0.3, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);
    
    // Test atoms from included files
    ASSERT_TRUE(parser.get_atom_properties("BENX", "CG", charge, mass));
    EXPECT_NEAR(charge, -0.115, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    ASSERT_TRUE(parser.get_atom_properties("PRPX", "C1", charge, mass));
    EXPECT_NEAR(charge, -0.27, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    ASSERT_TRUE(parser.get_atom_properties("SOL", "OW", charge, mass));
    EXPECT_NEAR(charge, -0.834, 1e-6);
    EXPECT_NEAR(mass, 15.9994, 1e-6);
}

TEST_F(TopParserTest, ParseTopologyWithMissingIncludes) {
    TopParser parser;
    std::string top_file = data_dir_ + "/test_missing_includes.top";
    
    // Should not fail even with missing includes
    ASSERT_TRUE(parser.parse_with_includes(top_file)) << "Parser failed with missing includes";
    
    // Should still be able to read atoms from main file
    double charge, mass;
    ASSERT_TRUE(parser.get_atom_properties("ALA", "N", charge, mass));
    EXPECT_NEAR(charge, -0.3, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);
}

TEST_F(TopParserTest, CompareWithAndWithoutIncludes) {
    TopParser parser_with_includes;
    TopParser parser_without_includes;
    std::string top_file = data_dir_ + "/test.top";
    
    ASSERT_TRUE(parser_with_includes.parse_with_includes(top_file)) << "Failed to parse with includes";
    ASSERT_TRUE(parser_without_includes.parse(top_file)) << "Failed to parse without includes";
    
    double charge1, mass1, charge2, mass2;
    
    // Main file atoms should be identical
    ASSERT_TRUE(parser_with_includes.get_atom_properties("ALA", "N", charge1, mass1));
    ASSERT_TRUE(parser_without_includes.get_atom_properties("ALA", "N", charge2, mass2));
    EXPECT_NEAR(charge1, charge2, 1e-6);
    EXPECT_NEAR(mass1, mass2, 1e-6);
    
    // Included atoms should only be available with includes enabled
    ASSERT_TRUE(parser_with_includes.get_atom_properties("SOL", "OW", charge1, mass1));
    ASSERT_FALSE(parser_without_includes.get_atom_properties("SOL", "OW", charge2, mass2));
}