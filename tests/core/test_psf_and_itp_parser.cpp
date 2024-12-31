// tests/core/test_psf_and_itp_parser.cpp

#include <gtest/gtest.h>
#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/io/itp_parser.hpp"
#include "config.hpp"
#include <filesystem>

using namespace pygcmc::core::io;

class PSFAndITPParserTest : public ::testing::Test {
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

TEST_F(PSFAndITPParserTest, ParsePSFAndITPFiles) {
    // First parse PSF file
    PSFParser psf_parser;
    std::string psf_file = data_dir_ + "/test_proa.psf";
    ASSERT_TRUE(psf_parser.parse(psf_file)) << "Failed to parse PSF file";

    // Then parse ITP files
    ITPParser itp_parser;
    std::string benx_file = data_dir_ + "/mols/benx.itp";
    std::string prpx_file = data_dir_ + "/mols/prpx.itp";
    std::string sol_file = data_dir_ + "/mols/sol.itp";

    ASSERT_TRUE(itp_parser.parse(benx_file)) << "Failed to parse BENX ITP file";
    ASSERT_TRUE(itp_parser.parse(prpx_file)) << "Failed to parse PRPX ITP file";
    ASSERT_TRUE(itp_parser.parse(sol_file)) << "Failed to parse SOL ITP file";

    // Test getting atom properties from PSF file
    double charge, mass;
    // Test ALA residue atoms (N-terminal, residue 7)
    ASSERT_TRUE(psf_parser.get_atom_properties("ALA", 7, "N", charge, mass));
    EXPECT_NEAR(charge, -0.3, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);

    // Test VAL residue atoms (residue 8)
    ASSERT_TRUE(psf_parser.get_atom_properties("VAL", 8, "CA", charge, mass));
    EXPECT_NEAR(charge, 0.07, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    // Test PRO residue atoms (residue 9)
    ASSERT_TRUE(psf_parser.get_atom_properties("PRO", 9, "N", charge, mass));
    EXPECT_NEAR(charge, -0.29, 1e-6);
    EXPECT_NEAR(mass, 14.007, 1e-6);

    // Test getting atom properties from ITP files
    // Test BENX residue atoms
    ASSERT_TRUE(itp_parser.get_atom_properties("BENX", "CG", charge, mass));
    EXPECT_NEAR(charge, -0.115, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    // Test PRPX residue atoms
    ASSERT_TRUE(itp_parser.get_atom_properties("PRPX", "C1", charge, mass));
    EXPECT_NEAR(charge, -0.27, 1e-6);
    EXPECT_NEAR(mass, 12.011, 1e-6);

    // Test SOL residue atoms
    ASSERT_TRUE(itp_parser.get_atom_properties("SOL", "OW", charge, mass));
    EXPECT_NEAR(charge, -0.834, 1e-6);
    EXPECT_NEAR(mass, 15.9994, 1e-6);
}

TEST_F(PSFAndITPParserTest, UpdatePDBAtoms) {
    // First parse PSF file
    PSFParser psf_parser;
    std::string psf_file = data_dir_ + "/test_proa.psf";
    ASSERT_TRUE(psf_parser.parse(psf_file)) << "Failed to parse PSF file";

    // Then parse ITP files
    ITPParser itp_parser;
    std::string benx_file = data_dir_ + "/mols/benx.itp";
    std::string prpx_file = data_dir_ + "/mols/prpx.itp";
    std::string sol_file = data_dir_ + "/mols/sol.itp";

    ASSERT_TRUE(itp_parser.parse(benx_file)) << "Failed to parse BENX ITP file";
    ASSERT_TRUE(itp_parser.parse(prpx_file)) << "Failed to parse PRPX ITP file";
    ASSERT_TRUE(itp_parser.parse(sol_file)) << "Failed to parse SOL ITP file";

    // Create test PDB atoms
    std::vector<PDBAtom> atoms;
    
    // Add atom from PSF
    PDBAtom ala_atom;
    ala_atom.residue = "ALA";
    ala_atom.name = "N";
    ala_atom.sequence = 1;
    atoms.push_back(ala_atom);

    // Add atom from ITP
    PDBAtom benx_atom;
    benx_atom.residue = "BENX";
    benx_atom.name = "CG";
    benx_atom.sequence = 1;
    atoms.push_back(benx_atom);

    // Update atoms from PSF
    int psf_updated = psf_parser.update_pdb_atoms(atoms);
    EXPECT_EQ(psf_updated, 1) << "Failed to update PDB atoms from PSF";
    EXPECT_NEAR(atoms[0].topo_charge, -0.3, 1e-6);
    EXPECT_NEAR(atoms[0].topo_mass, 14.007, 1e-6);

    // Update atoms from ITP
    int itp_updated = itp_parser.update_pdb_atoms(atoms);
    EXPECT_EQ(itp_updated, 1) << "Failed to update PDB atoms from ITP";
    EXPECT_NEAR(atoms[1].topo_charge, -0.115, 1e-6);
    EXPECT_NEAR(atoms[1].topo_mass, 12.011, 1e-6);
}

TEST_F(PSFAndITPParserTest, GetMissingTopologyInfo) {
    // First parse PSF file
    PSFParser psf_parser;
    std::string psf_file = data_dir_ + "/test_proa.psf";
    ASSERT_TRUE(psf_parser.parse(psf_file)) << "Failed to parse PSF file";

    // Then parse ITP files
    ITPParser itp_parser;
    std::string benx_file = data_dir_ + "/mols/benx.itp";
    std::string prpx_file = data_dir_ + "/mols/prpx.itp";
    std::string sol_file = data_dir_ + "/mols/sol.itp";

    ASSERT_TRUE(itp_parser.parse(benx_file)) << "Failed to parse BENX ITP file";
    ASSERT_TRUE(itp_parser.parse(prpx_file)) << "Failed to parse PRPX ITP file";
    ASSERT_TRUE(itp_parser.parse(sol_file)) << "Failed to parse SOL ITP file";

    // Create test PDB atoms with some missing topology
    std::vector<PDBAtom> atoms;
    
    // Add existing atoms
    PDBAtom ala_atom;
    ala_atom.residue = "ALA";
    ala_atom.name = "N";
    ala_atom.sequence = 1;
    atoms.push_back(ala_atom);

    PDBAtom benx_atom;
    benx_atom.residue = "BENX";
    benx_atom.name = "CG";
    benx_atom.sequence = 1;
    atoms.push_back(benx_atom);

    // Add non-existent atoms
    PDBAtom unknown_atom;
    unknown_atom.residue = "UNKNOWN";
    unknown_atom.name = "X";
    unknown_atom.sequence = 1;
    atoms.push_back(unknown_atom);

    // Update atoms from both parsers
    psf_parser.update_pdb_atoms(atoms);
    itp_parser.update_pdb_atoms(atoms);

    // Check missing topology info from PSF
    auto psf_missing = psf_parser.get_missing_topology_info(atoms);
    EXPECT_TRUE(psf_missing.find("BENX") != psf_missing.end()) << "BENX should be missing from PSF";
    EXPECT_TRUE(psf_missing.find("UNKNOWN") != psf_missing.end()) << "UNKNOWN should be missing from PSF";

    // Check missing topology info from ITP
    auto itp_missing = itp_parser.get_missing_topology_info(atoms);
    EXPECT_TRUE(itp_missing.find("ALA") != itp_missing.end()) << "ALA should be missing from ITP";
    EXPECT_TRUE(itp_missing.find("UNKNOWN") != itp_missing.end()) << "UNKNOWN should be missing from ITP";
}
