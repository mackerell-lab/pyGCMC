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
std::string g_water_ions_str = TEST_DATA_DIR "/toppar_water_ions.str";
std::string g_silcs_str = TEST_DATA_DIR "/silcs.str";

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
            } else if (arg.rfind("--water_ions_str=", 0) == 0) {
                g_water_ions_str = arg.substr(17);
            } else if (arg.rfind("--silcs_str=", 0) == 0) {
                g_silcs_str = arg.substr(12);
            }
        }

        if (g_pdb_file.empty() || g_top_file.empty() || g_cgenff_prm_file.empty() || 
            g_prot_prm_file.empty() || g_water_ions_str.empty() || g_silcs_str.empty()) {
            std::cerr << "Error: --pdb_file, --top_file, --cgenff_prm_file, --prot_prm_file, --water_ions_str, and --silcs_str arguments are required" << std::endl;
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
        ff_parser_water_ = std::make_unique<FFParser>();
        ff_parser_silcs_ = std::make_unique<FFParser>();
    }

    std::unique_ptr<PDBParser> pdb_parser_;
    std::unique_ptr<TopParser> top_parser_;
    std::unique_ptr<FFParser> ff_parser_;
    std::unique_ptr<FFParser> ff_parser_prot_;
    std::unique_ptr<FFParser> ff_parser_water_;
    std::unique_ptr<FFParser> ff_parser_silcs_;
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
    ASSERT_TRUE(ff_parser_water_->parse(g_water_ions_str))
        << "Failed to parse water and ions force field file";
    ASSERT_TRUE(ff_parser_silcs_->parse(g_silcs_str))
        << "Failed to parse SILCS force field file";

    // Debug: Print SILCS force field parameters
    const auto& silcs_params = ff_parser_silcs_->get_nonbonded_params();
    {
        auto it = silcs_params.find("LP");
        if (it != silcs_params.end()) {
            std::cout << "Found LP parameters: epsilon=" << it->second.epsilon 
                     << ", rmin=" << it->second.rmin << std::endl;
        } else {
            std::cout << "LP parameters not found!" << std::endl;
        }
        
        it = silcs_params.find("LQ");
        if (it != silcs_params.end()) {
            std::cout << "Found LQ parameters: epsilon=" << it->second.epsilon 
                     << ", rmin=" << it->second.rmin << std::endl;
        } else {
            std::cout << "LQ parameters not found!" << std::endl;
        }

        // Print NBFIX parameters
        const auto& nbfix_params = ff_parser_silcs_->get_nbfix_params();
        for (const auto& [key, value] : nbfix_params) {
            std::cout << "Found NBFIX parameters for " << key.first << "-" << key.second 
                     << ": epsilon=" << value.epsilon << ", rmin=" << value.rmin << std::endl;
        }
    }

    // Debug: Print water force field parameters
    const auto& water_params = ff_parser_water_->get_nonbonded_params();
    {
        auto it = water_params.find("OT");
        if (it != water_params.end()) {
            std::cout << "Found OT parameters: epsilon=" << it->second.epsilon 
                     << ", rmin=" << it->second.rmin << std::endl;
        } else {
            std::cout << "OT parameters not found!" << std::endl;
        }
        
        it = water_params.find("HT");
        if (it != water_params.end()) {
            std::cout << "Found HT parameters: epsilon=" << it->second.epsilon 
                     << ", rmin=" << it->second.rmin << std::endl;
        } else {
            std::cout << "HT parameters not found!" << std::endl;
        }
    }

    // 4. Update atoms with topology information
    int topo_updated = top_parser_->update_pdb_atoms(atoms);
    ASSERT_GT(topo_updated, 0) << "No atoms were updated with topology information";

    // Debug: Print topology types for water atoms
    for (const auto& atom : atoms) {
        if (atom.residue == "SOL") {
            std::cout << "Water atom: " << atom.name << ", topo_type=" << atom.topo_type 
                     << ", charge=" << atom.topo_charge << std::endl;
        }
    }

    // 5. Update atoms with force field parameters
    // Try each force field in sequence: protein, water/ions, then general
    int ff_updated_prot = ff_parser_prot_->update_pdb_atoms(atoms);
    int ff_updated_water = ff_parser_water_->update_pdb_atoms(atoms);
    int ff_updated_gen = ff_parser_->update_pdb_atoms(atoms);
    
    // Debug: Print updated force field parameters for water atoms
    for (const auto& atom : atoms) {
        if (atom.residue == "SOL") {
            std::cout << "Water atom after FF update: " << atom.name 
                     << ", epsilon=" << atom.forcefield_epsilon 
                     << ", rmin=" << atom.forcefield_rmin << std::endl;
        }
    }

    ASSERT_GT(ff_updated_prot + ff_updated_water + ff_updated_gen, 0) 
        << "No atoms were updated with force field parameters";

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
            EXPECT_NEAR(atom.forcefield_epsilon, -0.2, 1e-4);
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
            EXPECT_NEAR(atom.forcefield_epsilon, -0.032, 1e-4);
            EXPECT_NEAR(atom.forcefield_rmin, 4.0, 1e-4);
            break;
        }
    }
    EXPECT_TRUE(found_val_ca) << "Could not find VAL CA atom";

    // Test water oxygen atom
    bool found_water_o = false;
    for (const auto& atom : atoms) {
        if (atom.residue == "SOL" && atom.name == "OW") {
            found_water_o = true;
            
            // Check topology parameters
            EXPECT_NEAR(atom.topo_charge, -0.834, 1e-4);
            EXPECT_NEAR(atom.topo_mass, 15.9994, 1e-4);
            EXPECT_EQ(atom.topo_type, "OT");

            // Check force field parameters
            EXPECT_NEAR(atom.forcefield_epsilon, -0.1521, 1e-4);
            EXPECT_NEAR(atom.forcefield_rmin, 1.7682 * 2, 1e-4);
            break;
        }
    }
    EXPECT_TRUE(found_water_o) << "Could not find water oxygen atom";

    // Test water hydrogen atom
    bool found_water_h = false;
    for (const auto& atom : atoms) {
        if (atom.residue == "SOL" && atom.name == "HW1") {
            found_water_h = true;
            
            // Check topology parameters
            EXPECT_NEAR(atom.topo_charge, 0.417, 1e-4);
            EXPECT_NEAR(atom.topo_mass, 1.008, 1e-4);
            EXPECT_EQ(atom.topo_type, "HT");

            // Check force field parameters
            EXPECT_NEAR(atom.forcefield_epsilon, -0.046, 1e-4);
            EXPECT_NEAR(atom.forcefield_rmin, 0.2245 * 2, 1e-4);
            break;
        }
    }
    EXPECT_TRUE(found_water_h) << "Could not find water hydrogen atom";

    // Test total charge of water molecule
    double total_water_charge = 0.0;
    bool found_water = false;
    for (const auto& atom : atoms) {
        if (atom.residue == "SOL") {
            if (!found_water) {
                found_water = true;
                total_water_charge = 0.0;
            }
            total_water_charge += atom.topo_charge;
        } else if (found_water) {
            break;  // We've moved past the water molecule
        }
    }
    if (found_water) {
        EXPECT_NEAR(total_water_charge, 0.0, 1e-4) << "Water molecule should have neutral total charge";
    }

    // Print all force field parameters
    std::cout << "\n=== Force Field Parameters Summary ===\n";
    
    // Print global nonbonded parameters
    std::cout << "\nGlobal Nonbonded Parameters:\n";
    std::cout << "  cutnb  = " << ff_parser_->get_cutnb() << " Å\n";
    std::cout << "  ctofnb = " << ff_parser_->get_ctofnb() << " Å\n";
    std::cout << "  ctonnb = " << ff_parser_->get_ctonnb() << " Å\n";
    std::cout << "  eps    = " << ff_parser_->get_eps() << "\n";
    std::cout << "  e14fac = " << ff_parser_->get_e14fac() << "\n";
    std::cout << "  wmin   = " << ff_parser_->get_wmin() << "\n";

    // Print nonbonded parameters from each force field
    auto print_nonbonded = [](const std::string& name, const std::unique_ptr<FFParser>& parser) {
        const auto& params = parser->get_nonbonded_params();
        std::cout << "\nNonbonded Parameters from " << name << " (" 
                 << params.size() << " entries):\n";
        for (const auto& [type, param] : params) {
            std::cout << "  " << type << ": epsilon=" << param.epsilon 
                     << ", rmin=" << param.rmin << "\n";
        }
    };

    print_nonbonded("Protein FF", ff_parser_prot_);
    print_nonbonded("CGenFF", ff_parser_);
    print_nonbonded("Water/Ions", ff_parser_water_);
    print_nonbonded("SILCS", ff_parser_silcs_);

    // Print NBFIX parameters from each force field
    auto print_nbfix = [](const std::string& name, const std::unique_ptr<FFParser>& parser) {
        const auto& params = parser->get_nbfix_params();
        std::cout << "\nNBFIX Parameters from " << name << " (" 
                 << params.size() << " entries):\n";
        for (const auto& [types, param] : params) {
            std::cout << "  " << types.first << "-" << types.second 
                     << ": epsilon=" << param.epsilon 
                     << ", rmin=" << param.rmin << "\n";
        }
    };

    print_nbfix("Protein FF", ff_parser_prot_);
    print_nbfix("CGenFF", ff_parser_);
    print_nbfix("Water/Ions", ff_parser_water_);
    print_nbfix("SILCS", ff_parser_silcs_);

    std::cout << "\n=== End of Force Field Parameters ===\n";
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

