// tests/core/test_pdb_parser.cpp

#include "gtest/gtest.h"
#include "config.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <tuple>
#include <map>

using namespace pygcmc::core::io;

class PDBParserTest : public ::testing::Test {
protected:
    std::string test_pdb;
    std::string water_pdb;

    void SetUp() override {
        test_pdb = std::string(PDB_DATA_DIR) + "/test.pdb";
        water_pdb = std::string(PDB_DATA_DIR) + "/water.pdb";
    }

    void print_crystal_info(const std::string& filename, const std::optional<std::vector<double>>& cryst) {
        std::cout << "\nCrystal information for " << filename << ":\n";
        if (cryst) {
            const auto& box = cryst.value();
            std::cout << "Box dimensions found:\n";
            std::cout << "  a = " << box[0] << " Å\n";
            std::cout << "  b = " << box[1] << " Å\n";
            std::cout << "  c = " << box[2] << " Å\n";
            if (box.size() > 3) {
                std::cout << "  α = " << box[3] << "°\n";
                std::cout << "  β = " << box[4] << "°\n";
                std::cout << "  γ = " << box[5] << "°\n";
            }
        } else {
            std::cout << "No crystal information found.\n";
        }
        std::cout << std::endl;
    }

    void verify_atom_position(const PDBAtom& atom, double expected_x, 
                            double expected_y, double expected_z) {
        auto pos = atom.position();
        EXPECT_DOUBLE_EQ(pos[0], expected_x);
        EXPECT_DOUBLE_EQ(pos[1], expected_y);
        EXPECT_DOUBLE_EQ(pos[2], expected_z);
    }

    void verify_center_of_mass(const IOResidue& residue, double expected_x,
                             double expected_y, double expected_z) {
        auto com = residue.center_of_mass();
        EXPECT_NEAR(com[0], expected_x, 1e-10);
        EXPECT_NEAR(com[1], expected_y, 1e-10);
        EXPECT_NEAR(com[2], expected_z, 1e-10);
    }

    void print_residue_info(const IOResidue& residue) {
        auto com = residue.center_of_mass();
        std::cout << "Residue " << residue.sequence_number 
                  << " (" << residue.name << "):" << std::endl
                  << "  Atoms: " << residue.atom_count() << std::endl
                  << "  Center of mass: (" << com[0] << ", " 
                  << com[1] << ", " << com[2] << ")" << std::endl;
        
        for (const auto& atom : residue.atoms) {
            auto pos = atom.position();
            std::cout << "  Serial: " << atom.serial
                      << ", Name: " << atom.name
                      << ", Position: (" << pos[0] << ", " 
                      << pos[1] << ", " << pos[2] << ")"
                      << ", Element: " << atom.element
                      << ", Type: " << atom.type
                      << std::endl;
        }
    }
};

class PDBParserParamTest : public ::testing::TestWithParam<std::tuple<std::string, bool>> {
protected:
    std::string pdb_file;
    bool is_water;

    void SetUp() override {
        std::tie(pdb_file, is_water) = GetParam();
        std::cout << "Using PDB file: " << pdb_file << ", is_water: " << is_water << std::endl;
    }

    void verify_atom_position(const PDBAtom& atom, double expected_x, 
                            double expected_y, double expected_z) {
        auto pos = atom.position();
        EXPECT_DOUBLE_EQ(pos[0], expected_x);
        EXPECT_DOUBLE_EQ(pos[1], expected_y);
        EXPECT_DOUBLE_EQ(pos[2], expected_z);
    }

    void verify_center_of_mass(const IOResidue& residue, double expected_x,
                             double expected_y, double expected_z) {
        auto com = residue.center_of_mass();
        EXPECT_NEAR(com[0], expected_x, 1e-10);
        EXPECT_NEAR(com[1], expected_y, 1e-10);
        EXPECT_NEAR(com[2], expected_z, 1e-10);
    }

    void print_residue_info(const IOResidue& residue) {
        auto com = residue.center_of_mass();
        std::cout << "Residue " << residue.sequence_number 
                  << " (" << residue.name << "):" << std::endl
                  << "  Atoms: " << residue.atom_count() << std::endl
                  << "  Center of mass: (" << com[0] << ", " 
                  << com[1] << ", " << com[2] << ")" << std::endl;
        
        for (const auto& atom : residue.atoms) {
            auto pos = atom.position();
            std::cout << "  Serial: " << atom.serial
                      << ", Name: " << atom.name
                      << ", Position: (" << pos[0] << ", " 
                      << pos[1] << ", " << pos[2] << ")"
                      << ", Element: " << atom.element
                      << ", Type: " << atom.type
                      << std::endl;
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    PDBFiles,
    PDBParserParamTest,
    ::testing::Values(
        std::make_tuple(std::string(PDB_DATA_DIR) + "/test.pdb", false),
        std::make_tuple(std::string(PDB_DATA_DIR) + "/water.pdb", true)
    )
);

TEST_P(PDBParserParamTest, ValidateResidueGeometry) {
    try {
        auto parsed = PDBParser::parse(pdb_file);
        
        if (is_water) {
            // Water molecule validation
            ASSERT_EQ(parsed.second.size(), 1);
            const auto& water = parsed.second[0];
            EXPECT_EQ(water.atom_count(), 3);
            
            // Verify water molecule geometry
            double expected_x = (0.0 + 0.957 - 0.24) / 3.0;
            double expected_y = (0.0 + 0.0 + 0.927) / 3.0;
            double expected_z = 0.0;
            verify_center_of_mass(water, expected_x, expected_y, expected_z);
            
            verify_atom_position(water.atoms[0], 0.0, 0.0, 0.0);      // O
            verify_atom_position(water.atoms[1], 0.957, 0.0, 0.0);    // H1
            verify_atom_position(water.atoms[2], -0.24, 0.927, 0.0);  // H2
        } else {
            // Regular PDB validation
            EXPECT_GT(parsed.second.size(), 1);
            for (const auto& residue : parsed.second) {
                EXPECT_TRUE(residue.is_valid());
            }
        }
    } catch (const ParserError& e) {
        FAIL() << "ParserError 异常: " << e.what();
    }
}

TEST_P(PDBParserParamTest, ValidateAtomPositions) {
    try {
        auto parsed = PDBParser::parse(pdb_file);
        for (const auto& residue : parsed.second) {
            for (const auto& atom : residue.atoms) {
                if (is_water) {
                    // Water atom positions
                    if (atom.name == "O") {
                        verify_atom_position(atom, 0.0, 0.0, 0.0);
                    } else if (atom.name == "H1") {
                        verify_atom_position(atom, 0.957, 0.0, 0.0);
                    } else if (atom.name == "H2") {
                        verify_atom_position(atom, -0.24, 0.927, 0.0);
                    }
                } else {
                    // Regular PDB atom positions
                    if (residue.name == "ALA" && residue.sequence_number == 7 && atom.name == "N") {
                        verify_atom_position(atom, 76.563, 93.118, 93.806);
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        FAIL() << "验证原子位置时发生异常: " << e.what();
    }
}

// Basic file handling tests
TEST_F(PDBParserTest, FileHandling) {
    EXPECT_THROW(PDBParser::parse("nonexistent.pdb"), ParserError);
    EXPECT_NO_THROW(PDBParser::parse(test_pdb));
    EXPECT_NO_THROW(PDBParser::parse(water_pdb));
}

// Crystal information tests
TEST_F(PDBParserTest, CrystalParameters) {
    // Test file with CRYST1 record
    auto parsed_test = PDBParser::parse(test_pdb);
    print_crystal_info("test.pdb", parsed_test.first);
    ASSERT_TRUE(parsed_test.first.has_value()) << "test.pdb should have crystal information";
    const auto& box = parsed_test.first.value();
    ASSERT_EQ(box.size(), 6) << "Box should have 6 parameters (a, b, c, alpha, beta, gamma)";
    EXPECT_DOUBLE_EQ(box[0], 127.022) << "Incorrect a parameter";
    EXPECT_DOUBLE_EQ(box[1], 133.419) << "Incorrect b parameter";
    EXPECT_DOUBLE_EQ(box[2], 132.854) << "Incorrect c parameter";
    EXPECT_DOUBLE_EQ(box[3], 90.0) << "Incorrect alpha angle";
    EXPECT_DOUBLE_EQ(box[4], 90.0) << "Incorrect beta angle";
    EXPECT_DOUBLE_EQ(box[5], 90.0) << "Incorrect gamma angle";

    // Test file with CRYST1 record (water.pdb)
    auto parsed_water = PDBParser::parse(water_pdb);
    print_crystal_info("water.pdb", parsed_water.first);
    ASSERT_TRUE(parsed_water.first.has_value()) << "water.pdb should have crystal information";
    const auto& water_box = parsed_water.first.value();
    ASSERT_EQ(water_box.size(), 6) << "Water box should have 6 parameters";
    EXPECT_DOUBLE_EQ(water_box[0], 10.0) << "Incorrect water box a parameter";
    EXPECT_DOUBLE_EQ(water_box[1], 10.0) << "Incorrect water box b parameter";
    EXPECT_DOUBLE_EQ(water_box[2], 10.0) << "Incorrect water box c parameter";
    EXPECT_DOUBLE_EQ(water_box[3], 90.0) << "Incorrect water box alpha angle";
    EXPECT_DOUBLE_EQ(water_box[4], 90.0) << "Incorrect water box beta angle";
    EXPECT_DOUBLE_EQ(water_box[5], 90.0) << "Incorrect water box gamma angle";

    // Test file without CRYST1 record (benx.pdb)
    std::string benx_pdb = std::string(PDB_DATA_DIR) + "/mols/benx.pdb";
    auto parsed_benx = PDBParser::parse(benx_pdb);
    print_crystal_info("benx.pdb", parsed_benx.first);
    EXPECT_FALSE(parsed_benx.first.has_value()) << "benx.pdb should not have crystal information";
    EXPECT_FALSE(parsed_benx.second.empty()) << "benx.pdb should still have residues";
}

// Residue parsing tests
TEST_F(PDBParserTest, ResidueStructure) {
    auto parsed = PDBParser::parse(test_pdb);
    const auto& residues = parsed.second;

    // Check first ALA residue (residue 7)
    const auto& ala = residues[0];
    EXPECT_EQ(ala.name, "ALA");
    EXPECT_EQ(ala.sequence_number, 7);
    EXPECT_EQ(ala.atom_count(), 12);
    
    // Check atom details for ALA
    const auto& n_atom = ala.atoms[0];
    EXPECT_EQ(n_atom.name, "N");
    EXPECT_EQ(n_atom.serial, 1);
    EXPECT_EQ(n_atom.element, "N");
    verify_atom_position(n_atom, 76.563, 93.118, 93.806);

    // Check last residue (should be water)
    const auto& last_water = residues.back();
    EXPECT_EQ(last_water.name, "SOL");
    EXPECT_EQ(last_water.sequence_number, 2910);
    EXPECT_EQ(last_water.atom_count(), 3);
}

// Water molecule specific tests
TEST_F(PDBParserTest, WaterMolecule) {
    auto parsed = PDBParser::parse(water_pdb);
    ASSERT_EQ(parsed.second.size(), 1);
    
    const auto& water = parsed.second[0];
    EXPECT_EQ(water.name, "HOH");
    EXPECT_EQ(water.sequence_number, 1);
    EXPECT_EQ(water.atom_count(), 3);
    
    // Check water geometry
    EXPECT_EQ(water.atoms[0].name, "O");
    EXPECT_EQ(water.atoms[1].name, "H1");
    EXPECT_EQ(water.atoms[2].name, "H2");
    
    verify_atom_position(water.atoms[0], 0.0, 0.0, 0.0);
    verify_atom_position(water.atoms[1], 0.957, 0.0, 0.0);
    verify_atom_position(water.atoms[2], -0.24, 0.927, 0.0);
}

// Chain ID and insertion code tests
TEST_F(PDBParserTest, ChainAndInsertionCodes) {
    auto parsed = PDBParser::parse(test_pdb);
    
    // Check chain continuity
    char current_chain = parsed.second[0].chain_id;
    for (const auto& residue : parsed.second) {
        if (residue.chain_id != current_chain) {
            // New chain found
            current_chain = residue.chain_id;
        }
        // Verify chain ID is valid
        EXPECT_TRUE(std::isalpha(residue.chain_id) || residue.chain_id == ' ');
    }
}

// Element derivation tests
TEST_F(PDBParserTest, ElementDerivation) {
    auto parsed = PDBParser::parse(test_pdb);
    
    std::map<std::string, std::string> expected_elements = {
        {"N", "N"}, {"H1", "H"}, {"CA", "C"}, {"O", "O"},
        {"CB", "C"}, {"CG", "C"}, {"CD", "C"}
    };
    
    for (const auto& residue : parsed.second) {
        for (const auto& atom : residue.atoms) {
            if (expected_elements.count(atom.name)) {
                EXPECT_EQ(atom.element, expected_elements[atom.name])
                    << "Atom name: " << atom.name;
            }
        }
    }
}

// Alternative location indicator tests
TEST_F(PDBParserTest, AlternativeLocations) {
    auto parsed = PDBParser::parse(test_pdb);
    
    for (const auto& residue : parsed.second) {
        for (const auto& atom : residue.atoms) {
            // Alt loc should be either space or A-Z
            EXPECT_TRUE(atom.alt_loc == ' ' || (atom.alt_loc >= 'A' && atom.alt_loc <= 'Z'))
                << "Invalid alt_loc: " << atom.alt_loc;
        }
    }
}
