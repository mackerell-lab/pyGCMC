// tests/core/test_psf_parser.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/io/parser_common.hpp"
#include <string>
#include <cmath>
#include <filesystem>
#include <fstream>

namespace pygcmc {
namespace core {
namespace io {
namespace test {

// Global variable for file path
std::string g_psf_file = TEST_DATA_DIR "/test_proa.psf";

// Parse command line arguments
class Environment : public ::testing::Environment {
public:
    ~Environment() override {}

    void SetUp() override {
        const auto& args = ::testing::internal::GetArgvs();
        
        for (const auto& arg : args) {
            if (arg.rfind("--psf_file=", 0) == 0) {
                g_psf_file = arg.substr(11);
            }
        }

        if (g_psf_file.empty()) {
            std::cerr << "Error: --psf_file argument is required" << std::endl;
            exit(1);
        }
    }

    void TearDown() override {}
};

TEST(PSFParserTest, DefaultTopologyValues) {
    // Test default values
    PDBAtom atom;
    EXPECT_TRUE(atom.topo_type.empty());
    EXPECT_TRUE(std::isnan(atom.topo_charge));
    EXPECT_TRUE(std::isnan(atom.topo_mass));
    EXPECT_FALSE(atom.has_topology_info());

    // Test partial topology information
    atom.topo_type = "CT1";
    EXPECT_FALSE(atom.has_topology_info());  // still false because charge and mass are NaN

    atom.topo_charge = -0.3;
    EXPECT_FALSE(atom.has_topology_info());  // still false because mass is NaN

    atom.topo_mass = 12.011;
    EXPECT_TRUE(atom.has_topology_info());  // now true because all fields are set

    // Test resetting to default
    atom.topo_type = "";
    atom.topo_charge = std::numeric_limits<double>::quiet_NaN();
    atom.topo_mass = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(atom.has_topology_info());
}

class TestPSFParser : public ::testing::Test {
protected:
    void SetUp() override {
        parser_ = std::make_unique<PSFParser>();
    }

    std::unique_ptr<PSFParser> parser_;
};

TEST_F(TestPSFParser, ParsePSFFile) {
    ASSERT_TRUE(parser_->parse(g_psf_file));
}

TEST_F(TestPSFParser, GetAtomProperties) {
    ASSERT_TRUE(parser_->parse(g_psf_file));

    double charge, mass;
    // Test ALA residue atoms (non N-terminal, residue 10)
    ASSERT_TRUE(parser_->get_atom_properties("ALA", 10, "N", charge, mass));
    EXPECT_NEAR(charge, -0.47, 1e-6);  // Non N-terminal ALA N atom
    EXPECT_NEAR(mass, 14.007, 1e-6);

    // Test VAL residue atoms (residue 8)
    ASSERT_TRUE(parser_->get_atom_properties("VAL", 8, "CA", charge, mass));
    EXPECT_NEAR(charge, 0.07, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    ASSERT_TRUE(parser_->get_atom_properties("VAL", 8, "CG1", charge, mass));
    EXPECT_NEAR(charge, -0.27, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    // Test PRO residue atoms (residue 9)
    ASSERT_TRUE(parser_->get_atom_properties("PRO", 9, "N", charge, mass));
    EXPECT_NEAR(charge, -0.29, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);

    ASSERT_TRUE(parser_->get_atom_properties("PRO", 9, "CD", charge, mass));
    EXPECT_NEAR(charge, 0.0, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    // Test ASN residue atoms (residue 12)
    ASSERT_TRUE(parser_->get_atom_properties("ASN", 12, "CG", charge, mass));
    EXPECT_NEAR(charge, 0.55, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    ASSERT_TRUE(parser_->get_atom_properties("ASN", 12, "OD1", charge, mass));
    EXPECT_NEAR(charge, -0.55, 1e-6);
    EXPECT_NEAR(mass, 15.9994, 1e-6);

    // Test GLN residue atoms (residue 13)
    ASSERT_TRUE(parser_->get_atom_properties("GLN", 13, "NE2", charge, mass));
    EXPECT_NEAR(charge, -0.62, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);

    ASSERT_TRUE(parser_->get_atom_properties("GLN", 13, "HE21", charge, mass));
    EXPECT_NEAR(charge, 0.32, 1e-6);
    EXPECT_NEAR(mass, 1.008, 1e-6);

    // Test ALA residue atom properties
    ASSERT_TRUE(parser_->get_atom_properties("ALA", 7, "HT1", charge, mass));
    EXPECT_NEAR(charge, 0.33, 1e-6);
    EXPECT_NEAR(mass, 1.008, 1e-6);

    // Test non-existent atoms and residues
    ASSERT_FALSE(parser_->get_atom_properties("XXX", 1, "XXX", charge, mass));
    ASSERT_FALSE(parser_->get_atom_properties("BEN", 1, "C1", charge, mass));
    ASSERT_FALSE(parser_->get_atom_properties("SOL", "OW", charge, mass));
}

TEST_F(TestPSFParser, UpdatePDBAtoms) {
    ASSERT_TRUE(parser_->parse(g_psf_file));

    std::vector<PDBAtom> pdb_atoms;
    
    // Create test PDB atoms
    // ALA residue atoms (N-terminal)
    PDBAtom ala_n(1, "N", "ALA", 7);
    pdb_atoms.push_back(ala_n);

    PDBAtom ala_ht1(2, "HT1", "ALA", 7);
    pdb_atoms.push_back(ala_ht1);

    PDBAtom ala_ca(3, "CA", "ALA", 7);
    pdb_atoms.push_back(ala_ca);

    PDBAtom ala_cb(4, "CB", "ALA", 7);
    pdb_atoms.push_back(ala_cb);

    // VAL residue atoms
    PDBAtom val_n(5, "N", "VAL", 8);
    pdb_atoms.push_back(val_n);

    PDBAtom val_ca(6, "CA", "VAL", 8);
    pdb_atoms.push_back(val_ca);

    PDBAtom val_cg1(7, "CG1", "VAL", 8);
    pdb_atoms.push_back(val_cg1);

    // PRO residue atoms
    PDBAtom pro_n(8, "N", "PRO", 9);
    pdb_atoms.push_back(pro_n);

    PDBAtom pro_cd(9, "CD", "PRO", 9);
    pdb_atoms.push_back(pro_cd);

    // ASN residue atoms
    PDBAtom asn_cg(10, "CG", "ASN", 10);
    pdb_atoms.push_back(asn_cg);

    PDBAtom asn_od1(11, "OD1", "ASN", 10);
    pdb_atoms.push_back(asn_od1);

    // GLN residue atoms
    PDBAtom gln_ne2(12, "NE2", "GLN", 11);
    pdb_atoms.push_back(gln_ne2);

    PDBAtom gln_he21(13, "HE21", "GLN", 11);
    pdb_atoms.push_back(gln_he21);

    // Update atoms with PSF topology
    int updated = parser_->update_pdb_atoms(pdb_atoms);
    EXPECT_EQ(updated, 13);

    // Verify ALA residue atoms (N-terminal)
    EXPECT_NEAR(pdb_atoms[0].topo_charge, -0.3, 1e-6);
    EXPECT_NEAR(pdb_atoms[0].topo_mass, 14.007, 1e-6);
    EXPECT_EQ(pdb_atoms[0].topo_type, "NH3");

    EXPECT_NEAR(pdb_atoms[1].topo_charge, 0.33, 1e-6);
    EXPECT_NEAR(pdb_atoms[1].topo_mass, 1.008, 1e-6);
    EXPECT_EQ(pdb_atoms[1].topo_type, "HC");

    EXPECT_NEAR(pdb_atoms[2].topo_charge, 0.21, 1e-6);
    EXPECT_NEAR(pdb_atoms[2].topo_mass, 12.011, 1e-6);
    EXPECT_EQ(pdb_atoms[2].topo_type, "CT1");

    EXPECT_NEAR(pdb_atoms[3].topo_charge, -0.27, 1e-6);
    EXPECT_NEAR(pdb_atoms[3].topo_mass, 12.011, 1e-6);
    EXPECT_EQ(pdb_atoms[3].topo_type, "CT3");

    // Verify VAL residue atoms
    EXPECT_NEAR(pdb_atoms[4].topo_charge, -0.47, 1e-6);
    EXPECT_NEAR(pdb_atoms[4].topo_mass, 14.007, 1e-6);
    EXPECT_EQ(pdb_atoms[4].topo_type, "NH1");

    EXPECT_NEAR(pdb_atoms[5].topo_charge, 0.07, 1e-6);
    EXPECT_NEAR(pdb_atoms[5].topo_mass, 12.011, 1e-6);
    EXPECT_EQ(pdb_atoms[5].topo_type, "CT1");

    EXPECT_NEAR(pdb_atoms[6].topo_charge, -0.27, 1e-6);
    EXPECT_NEAR(pdb_atoms[6].topo_mass, 12.011, 1e-6);
    EXPECT_EQ(pdb_atoms[6].topo_type, "CT3");

    // Verify PRO residue atoms
    EXPECT_NEAR(pdb_atoms[7].topo_charge, -0.29, 1e-6);
    EXPECT_NEAR(pdb_atoms[7].topo_mass, 14.007, 1e-6);
    EXPECT_EQ(pdb_atoms[7].topo_type, "N");

    EXPECT_NEAR(pdb_atoms[8].topo_charge, 0.0, 1e-6);
    EXPECT_NEAR(pdb_atoms[8].topo_mass, 12.011, 1e-6);
    EXPECT_EQ(pdb_atoms[8].topo_type, "CP3");

    // Verify ASN residue atoms
    EXPECT_NEAR(pdb_atoms[9].topo_charge, 0.55, 1e-6);
    EXPECT_NEAR(pdb_atoms[9].topo_mass, 12.011, 1e-6);
    EXPECT_EQ(pdb_atoms[9].topo_type, "CC");

    EXPECT_NEAR(pdb_atoms[10].topo_charge, -0.55, 1e-6);
    EXPECT_NEAR(pdb_atoms[10].topo_mass, 15.9994, 1e-6);
    EXPECT_EQ(pdb_atoms[10].topo_type, "O");

    // Verify GLN residue atoms
    EXPECT_NEAR(pdb_atoms[11].topo_charge, -0.62, 1e-6);
    EXPECT_NEAR(pdb_atoms[11].topo_mass, 14.007, 1e-6);
    EXPECT_EQ(pdb_atoms[11].topo_type, "NH2");

    EXPECT_NEAR(pdb_atoms[12].topo_charge, 0.32, 1e-6);
    EXPECT_NEAR(pdb_atoms[12].topo_mass, 1.008, 1e-6);
    EXPECT_EQ(pdb_atoms[12].topo_type, "H");
}

TEST_F(TestPSFParser, GetMissingTopologyInfo) {
    ASSERT_TRUE(parser_->parse(g_psf_file));

    std::vector<PDBAtom> pdb_atoms;
    
    // Add atoms with topology
    PDBAtom ala_n(1, "N", "ALA", 7);
    pdb_atoms.push_back(ala_n);

    PDBAtom val_ca(2, "CA", "VAL", 8);
    pdb_atoms.push_back(val_ca);

    // Add BEN (benzene) atoms without topology
    PDBAtom ben_cd1(3, "CD1", "BEN", 1);
    pdb_atoms.push_back(ben_cd1);

    PDBAtom ben_cd2(4, "CD2", "BEN", 1);
    pdb_atoms.push_back(ben_cd2);

    PDBAtom ben_ce1(5, "CE1", "BEN", 1);
    pdb_atoms.push_back(ben_ce1);

    PDBAtom ben_ce2(6, "CE2", "BEN", 1);
    pdb_atoms.push_back(ben_ce2);

    PDBAtom ben_cg(7, "CG", "BEN", 1);
    pdb_atoms.push_back(ben_cg);

    PDBAtom ben_cz(8, "CZ", "BEN", 1);
    pdb_atoms.push_back(ben_cz);

    // Add PRP (propane) atoms without topology
    PDBAtom prp_c1(9, "C1", "PRP", 2);
    pdb_atoms.push_back(prp_c1);

    PDBAtom prp_c2(10, "C2", "PRP", 2);
    pdb_atoms.push_back(prp_c2);

    PDBAtom prp_c3(11, "C3", "PRP", 2);
    pdb_atoms.push_back(prp_c3);

    // Add SOL (water) atoms without topology
    PDBAtom sol_ow(12, "OW", "SOL", 3);
    pdb_atoms.push_back(sol_ow);

    PDBAtom sol_hw1(13, "HW1", "SOL", 3);
    pdb_atoms.push_back(sol_hw1);

    PDBAtom sol_hw2(14, "HW2", "SOL", 3);
    pdb_atoms.push_back(sol_hw2);

    auto missing_info = parser_->get_missing_topology_info(pdb_atoms);

    // Verify expected missing residues
    EXPECT_TRUE(missing_info.find("BEN") != missing_info.end());
    EXPECT_TRUE(missing_info.find("PRP") != missing_info.end());
    EXPECT_TRUE(missing_info.find("SOL") != missing_info.end());

    // Verify that known residues have complete topology info
    EXPECT_TRUE(missing_info.find("ALA") == missing_info.end());
    EXPECT_TRUE(missing_info.find("VAL") == missing_info.end());
    EXPECT_TRUE(missing_info.find("PRO") == missing_info.end());
    EXPECT_TRUE(missing_info.find("ASN") == missing_info.end());
    EXPECT_TRUE(missing_info.find("GLN") == missing_info.end());

    // Verify BEN missing atoms
    if (missing_info.find("BEN") != missing_info.end()) {
        const auto& ben_atoms = missing_info.at("BEN");
        EXPECT_TRUE(ben_atoms.find("CD1") != ben_atoms.end());
        EXPECT_TRUE(ben_atoms.find("CD2") != ben_atoms.end());
        EXPECT_TRUE(ben_atoms.find("CE1") != ben_atoms.end());
        EXPECT_TRUE(ben_atoms.find("CE2") != ben_atoms.end());
        EXPECT_TRUE(ben_atoms.find("CG") != ben_atoms.end());
        EXPECT_TRUE(ben_atoms.find("CZ") != ben_atoms.end());
    }

    // Verify PRP missing atoms
    if (missing_info.find("PRP") != missing_info.end()) {
        const auto& prp_atoms = missing_info.at("PRP");
        EXPECT_TRUE(prp_atoms.find("C1") != prp_atoms.end());
        EXPECT_TRUE(prp_atoms.find("C2") != prp_atoms.end());
        EXPECT_TRUE(prp_atoms.find("C3") != prp_atoms.end());
    }

    // Verify SOL missing atoms
    if (missing_info.find("SOL") != missing_info.end()) {
        const auto& sol_atoms = missing_info.at("SOL");
        EXPECT_TRUE(sol_atoms.find("OW") != sol_atoms.end());
        EXPECT_TRUE(sol_atoms.find("HW1") != sol_atoms.end());
        EXPECT_TRUE(sol_atoms.find("HW2") != sol_atoms.end());
    }

    // Print missing topology information for debugging
    for (const auto& [residue, atoms] : missing_info) {
        std::cout << "Missing topology for " << residue << " atoms:";
        for (const auto& atom : atoms) {
            std::cout << " " << atom;
        }
        std::cout << std::endl;
    }
}

TEST_F(TestPSFParser, ParsePSFFileWithRoughMode) {
    // Test rough parsing mode
    ASSERT_TRUE(parser_->parse(g_psf_file, PSFParsingMode::Rough));

    // Verify that essential fields are parsed
    double charge, mass;
    ASSERT_TRUE(parser_->get_atom_properties("ALA", "N", charge, mass));
    // In rough mode, charge and mass should be 0
    EXPECT_NEAR(charge, 0.0, 1e-6);
    EXPECT_NEAR(mass, 0.0, 1e-6);
}

TEST_F(TestPSFParser, ParseMultiplePSFFiles) {
    std::vector<std::string> psf_files = {
        g_psf_file,  // Main PSF file
        TEST_DATA_DIR "/mols/sol.psf"  // Water PSF file
    };

    ASSERT_TRUE(parser_->parse_files(psf_files));

    // Test that atoms from both files are present
    double charge, mass;
    
    // Check protein atoms (from main PSF)
    ASSERT_TRUE(parser_->get_atom_properties("ALA", 7, "N", charge, mass));
    EXPECT_NEAR(charge, -0.3, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);

    // Check water atoms (from sol.psf)
    ASSERT_TRUE(parser_->get_atom_properties("SOL", "OW", charge, mass));
    EXPECT_NEAR(charge, -0.834, 1e-6);
    EXPECT_NEAR(mass, 15.9994, 1e-6);
}

TEST_F(TestPSFParser, HandleMalformedPSFFile) {
    // Create a temporary malformed PSF file
    std::string temp_psf = "temp_malformed.psf";
    {
        std::ofstream file(temp_psf);
        file << "PSF\n\n";
        file << "       1 !NATOM\n";
        file << "       1 MAIN 1    ALA  N    NH3   -0.30\n";  // Missing mass field
        file << "\n";
    }

    // Test that exact parsing fails
    EXPECT_THROW(parser_->parse(temp_psf, PSFParsingMode::Exact), std::runtime_error);

    // Test that rough parsing succeeds
    ASSERT_TRUE(parser_->parse(temp_psf, PSFParsingMode::Rough));

    // Verify the atom was parsed with default values
    double charge, mass;
    ASSERT_TRUE(parser_->get_atom_properties("ALA", "N", charge, mass));
    EXPECT_NEAR(charge, 0.0, 1e-6);
    EXPECT_NEAR(mass, 0.0, 1e-6);

    // Clean up
    std::filesystem::remove(temp_psf);
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