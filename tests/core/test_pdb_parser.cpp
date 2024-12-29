// tests/core/test_pdb_parser.cpp

#include "gtest/gtest.h"
#include "config.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <tuple>

using namespace pygcmc::core::io;

class PDBParserTest : public ::testing::Test {
protected:
    std::string pdb_file;

    void SetUp() override {
        // Get the PDB file path from command line arguments
        const auto args = ::testing::internal::GetArgvs();
        
        for (const auto& arg : args) {
            if (arg.find("--pdb_file=") == 0) {
                pdb_file = arg.substr(11);
                break;
            }
        }
        
        if (pdb_file.empty()) {
            throw std::runtime_error("PDB file path not provided. Use --pdb_file=<path>");
        }
        
        // Print the file path for debugging
        std::cout << "Using PDB file: " << pdb_file << std::endl;
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
