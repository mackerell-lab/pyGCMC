// tests/core/test_pdb_top_ff_parser.cpp

#include <gtest/gtest.h>
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/ff_parser.hpp"
#include "pygcmc/core/io/parser_common.hpp"
#include <string>
#include <cmath>
#include <memory>
#include <filesystem>

namespace pygcmc {
namespace core {
namespace io {
namespace test {

// Global variables for file paths
std::string g_pdb_file = TEST_DATA_DIR "/test.pdb";
std::string g_top_file = TEST_DATA_DIR "/test.top";
std::string g_cgenff_prm_file = TEST_DATA_DIR "/par_all36_cgenff.prm";
std::string g_prot_prm_file = TEST_DATA_DIR "/par_all36m_prot.prm";

// Parse command line arguments
class Environment : public ::testing::Environment {
public:
    ~Environment() override {}

    void SetUp() override {
        const auto& args = ::testing::internal::GetArgvs();
        
        for (const auto& arg : args) {
            if (arg.rfind("--pdb_file=", 0) == 0) {
                g_pdb_file = arg.substr(11);
            } else if (arg.rfind("--top_file=", 0) == 0) {
                g_top_file = arg.substr(11);
            } else if (arg.rfind("--cgenff_prm_file=", 0) == 0) {
                g_cgenff_prm_file = arg.substr(18);
            } else if (arg.rfind("--prot_prm_file=", 0) == 0) {
                g_prot_prm_file = arg.substr(16);
            }
        }

        if (g_pdb_file.empty() || g_top_file.empty() || g_cgenff_prm_file.empty() || g_prot_prm_file.empty()) {
            std::cerr << "Error: --pdb_file, --top_file, --cgenff_prm_file, and --prot_prm_file arguments are required" << std::endl;
            exit(1);
        }
    }

    void TearDown() override {}
};

class PDBTopFFParserTest : public ::testing::Test {
protected:
    void SetUp() override {
        pdb_parser_ = std::make_unique<PDBParser>();
        top_parser_ = std::make_unique<TopParser>();
        ff_parser_ = std::make_unique<FFParser>();
        ff_parser_prot_ = std::make_unique<FFParser>();
    }

    std::unique_ptr<PDBParser> pdb_parser_;
    std::unique_ptr<TopParser> top_parser_;
    std::unique_ptr<FFParser> ff_parser_;
    std::unique_ptr<FFParser> ff_parser_prot_;
};

TEST_F(PDBTopFFParserTest, CombinePDBTopFF) {
    // 1. Parse PDB file
    auto [coords, residues] = PDBParser::parse(g_pdb_file);
    ASSERT_FALSE(residues.empty()) << "PDB file contains no residues";

    // Extract atoms from residues
    std::vector<PDBAtom> atoms;
    for (const auto& residue : residues) {
        atoms.insert(atoms.end(), residue.atoms.begin(), residue.atoms.end());
    }

    // 2. Parse topology file with includes
    ASSERT_TRUE(top_parser_->parse_with_includes(g_top_file)) 
        << "Failed to parse topology file";

    // 3. Parse force field files
    ASSERT_TRUE(ff_parser_prot_->parse(g_prot_prm_file))
        << "Failed to parse protein force field file";
    ASSERT_TRUE(ff_parser_->parse(g_cgenff_prm_file))
        << "Failed to parse general force field file";

    // 4. Update atoms with topology information
    int topo_updated = top_parser_->update_pdb_atoms(atoms);
    ASSERT_GT(topo_updated, 0) << "No atoms were updated with topology information";

    // 5. Update atoms with force field parameters
    // First try protein force field, then general force field for any remaining atoms
    int ff_updated_prot = ff_parser_prot_->update_pdb_atoms(atoms);
    int ff_updated_gen = ff_parser_->update_pdb_atoms(atoms);
    ASSERT_GT(ff_updated_prot + ff_updated_gen, 0) << "No atoms were updated with force field parameters";

    // 6. Verify specific atoms
    // Test ALA N atom (residue 7)
    bool found_ala_n = false;
    for (const auto& atom : atoms) {
        if (atom.residue == "ALA" && atom.sequence == 7 && atom.name == "N") {
            found_ala_n = true;
            
            // Check topology parameters
            EXPECT_NEAR(atom.topo_charge, -0.3, 1e-4);
            EXPECT_NEAR(atom.topo_mass, 14.007, 1e-4);
            EXPECT_EQ(atom.topo_type, "NH3");

            // Check force field parameters - Rmin doubled from 1.85 to 3.70
            EXPECT_NEAR(atom.forcefield_epsilon, 0.2, 1e-4);
            EXPECT_NEAR(atom.forcefield_rmin, 3.70, 1e-4);
            break;
        }
    }
    EXPECT_TRUE(found_ala_n) << "Could not find ALA N atom";

    // Test VAL CA atom (residue 8)
    bool found_val_ca = false;
    for (const auto& atom : atoms) {
        if (atom.residue == "VAL" && atom.sequence == 8 && atom.name == "CA") {
            found_val_ca = true;
            
            // Check topology parameters
            EXPECT_NEAR(atom.topo_charge, 0.07, 1e-4);
            EXPECT_NEAR(atom.topo_mass, 12.011, 1e-4);
            EXPECT_EQ(atom.topo_type, "CT1");

            // Check force field parameters - Rmin doubled from 2.0 to 4.0
            EXPECT_NEAR(atom.forcefield_epsilon, 0.032, 1e-4);
            EXPECT_NEAR(atom.forcefield_rmin, 4.0, 1e-4);
            break;
        }
    }
    EXPECT_TRUE(found_val_ca) << "Could not find VAL CA atom";
}

} // namespace test
} // namespace io
} // namespace core
} // namespace pygcmc

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new pygcmc::core::io::test::Environment);
    return RUN_ALL_TESTS();
}

