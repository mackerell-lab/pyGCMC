// tests/core/test_psf_and_itp_parser.cpp

#include <gtest/gtest.h>
#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/io/itp_parser.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
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

TEST_F(PSFAndITPParserTest, CompareSolPSFAndITP) {
    // Parse SOL ITP file
    ITPParser itp_parser;
    std::string sol_itp_file = data_dir_ + "/mols/sol.itp";
    ASSERT_TRUE(itp_parser.parse(sol_itp_file)) << "Failed to parse SOL ITP file";

    // Parse SOL PSF file
    PSFParser psf_parser;
    std::string sol_psf_file = data_dir_ + "/mols/sol.psf";
    ASSERT_TRUE(psf_parser.parse(sol_psf_file)) << "Failed to parse SOL PSF file";

    // Create multiple water molecules with different sequence numbers
    std::vector<PDBAtom> water_atoms;
    
    // First water molecule (SOL1)
    PDBAtom ow1, hw11, hw12;
    ow1.residue = "SOL"; ow1.name = "OW"; ow1.sequence = 1;
    hw11.residue = "SOL"; hw11.name = "HW1"; hw11.sequence = 1;
    hw12.residue = "SOL"; hw12.name = "HW2"; hw12.sequence = 1;
    water_atoms.push_back(ow1);
    water_atoms.push_back(hw11);
    water_atoms.push_back(hw12);

    // Second water molecule (SOL2) - different sequence number
    PDBAtom ow2, hw21, hw22;
    ow2.residue = "SOL"; ow2.name = "OW"; ow2.sequence = 2;
    hw21.residue = "SOL"; hw21.name = "HW1"; hw21.sequence = 2;
    hw22.residue = "SOL"; hw22.name = "HW2"; hw22.sequence = 2;
    water_atoms.push_back(ow2);
    water_atoms.push_back(hw21);
    water_atoms.push_back(hw22);

    // Third water molecule (SOL3) - yet another sequence number
    PDBAtom ow3, hw31, hw32;
    ow3.residue = "SOL"; ow3.name = "OW"; ow3.sequence = 3;
    hw31.residue = "SOL"; hw31.name = "HW1"; hw31.sequence = 3;
    hw32.residue = "SOL"; hw32.name = "HW2"; hw32.sequence = 3;
    water_atoms.push_back(ow3);
    water_atoms.push_back(hw31);
    water_atoms.push_back(hw32);

    // First update using ITP parser
    std::vector<PDBAtom> itp_atoms = water_atoms;
    int itp_updated = itp_parser.update_pdb_atoms(itp_atoms);
    EXPECT_EQ(itp_updated, 9) << "Failed to update all water atoms from ITP";

    // Then update using PSF parser
    std::vector<PDBAtom> psf_atoms = water_atoms;
    int psf_updated = psf_parser.update_pdb_atoms(psf_atoms);
    EXPECT_EQ(psf_updated, 9) << "Failed to update all water atoms from PSF";

    // Compare results for each atom
    for (size_t i = 0; i < water_atoms.size(); ++i) {
        SCOPED_TRACE("Comparing atom " + std::to_string(i) + 
                     " (" + water_atoms[i].residue + " " + 
                     water_atoms[i].name + " " + 
                     std::to_string(water_atoms[i].sequence) + ")");
        
        // Compare topology type
        EXPECT_EQ(itp_atoms[i].topo_type, psf_atoms[i].topo_type)
            << "Topology type mismatch";
        
        // Compare charge (with small tolerance for floating point comparison)
        EXPECT_NEAR(itp_atoms[i].topo_charge, psf_atoms[i].topo_charge, 1e-6)
            << "Charge mismatch";
        
        // Compare mass (with small tolerance for floating point comparison)
        EXPECT_NEAR(itp_atoms[i].topo_mass, psf_atoms[i].topo_mass, 1e-6)
            << "Mass mismatch";

        // Verify expected values for water model
        if (water_atoms[i].name == "OW") {
            EXPECT_NEAR(psf_atoms[i].topo_charge, -0.834, 1e-6) << "Wrong charge for OW";
            EXPECT_NEAR(psf_atoms[i].topo_mass, 15.9994, 1e-6) << "Wrong mass for OW";
        } else if (water_atoms[i].name == "HW1" || water_atoms[i].name == "HW2") {
            EXPECT_NEAR(psf_atoms[i].topo_charge, 0.417, 1e-6) << "Wrong charge for HW";
            EXPECT_NEAR(psf_atoms[i].topo_mass, 1.008, 1e-6) << "Wrong mass for HW";
        }
    }
}

// TEST_F(PSFAndITPParserTest, CompareCompleteStructures) {
//     // Parse PSF files
//     PSFParser psf_parser;
//     std::string protein_psf_file = data_dir_ + "/test_proa.psf";
//     std::string sol_psf_file = data_dir_ + "/mols/sol.psf";
//     ASSERT_TRUE(psf_parser.parse(protein_psf_file)) << "Failed to parse protein PSF file";
//     ASSERT_TRUE(psf_parser.parse(sol_psf_file)) << "Failed to read SOL PSF file";

//     // Parse ITP files
//     ITPParser itp_parser;
//     std::string top_file = data_dir_ + "/test.top";
//     std::string sol_itp = data_dir_ + "/mols/sol.itp";
//     ASSERT_TRUE(itp_parser.parse(top_file)) << "Failed to parse TOP file";
//     ASSERT_TRUE(itp_parser.parse(sol_itp)) << "Failed to parse SOL ITP file";

//     // Create test atoms from PDB file
//     PDBParser pdb_parser;
//     std::string pdb_file = data_dir_ + "/test.pdb";
//     auto [box_info, residues] = pdb_parser.parse(pdb_file);
    
//     // Convert residues to a flat vector of atoms
//     std::vector<PDBAtom> pdb_atoms;
//     for (const auto& residue : residues) {
//         pdb_atoms.insert(pdb_atoms.end(), residue.atoms.begin(), residue.atoms.end());
//     }

//     // Create two copies of the atoms for PSF and ITP updates
//     std::vector<PDBAtom> psf_atoms = pdb_atoms;
//     std::vector<PDBAtom> itp_atoms = pdb_atoms;

//     // Update atoms using both parsers
//     int psf_updated = psf_parser.update_pdb_atoms(psf_atoms);
//     int itp_updated = itp_parser.update_pdb_atoms(itp_atoms);

//     // Compare number of updated atoms
//     EXPECT_EQ(psf_updated, itp_updated) << "Different number of atoms updated";

//     // Compare each atom's properties
//     ASSERT_EQ(psf_atoms.size(), itp_atoms.size()) << "Number of atoms mismatch";
//     for (size_t i = 0; i < psf_atoms.size(); ++i) {
//         SCOPED_TRACE("Comparing atom " + std::to_string(i) + 
//                      " (" + psf_atoms[i].residue + " " + 
//                      psf_atoms[i].name + " " + 
//                      std::to_string(psf_atoms[i].sequence) + ")");

//         // Compare basic properties
//         EXPECT_EQ(psf_atoms[i].residue, itp_atoms[i].residue) << "Residue mismatch";
//         EXPECT_EQ(psf_atoms[i].name, itp_atoms[i].name) << "Atom name mismatch";
//         EXPECT_EQ(psf_atoms[i].sequence, itp_atoms[i].sequence) << "Sequence number mismatch";

//         // Compare topology information
//         EXPECT_EQ(psf_atoms[i].topo_type, itp_atoms[i].topo_type) << "Topology type mismatch";
//         EXPECT_NEAR(psf_atoms[i].topo_charge, itp_atoms[i].topo_charge, 1e-6) << "Charge mismatch";
//         EXPECT_NEAR(psf_atoms[i].topo_mass, itp_atoms[i].topo_mass, 1e-6) << "Mass mismatch";

//         // Compare coordinates
//         EXPECT_NEAR(psf_atoms[i].x, itp_atoms[i].x, 1e-6) << "X coordinate mismatch";
//         EXPECT_NEAR(psf_atoms[i].y, itp_atoms[i].y, 1e-6) << "Y coordinate mismatch";
//         EXPECT_NEAR(psf_atoms[i].z, itp_atoms[i].z, 1e-6) << "Z coordinate mismatch";
//     }

//     // Compare box parameters
//     if (box_info.has_value()) {
//         const auto& box = box_info.value();
//         ASSERT_EQ(box.size(), 6) << "Box parameters should have 6 values";
        
//         // Verify box parameters are reasonable
//         EXPECT_GT(box[0], 0) << "Box length a should be positive";
//         EXPECT_GT(box[1], 0) << "Box length b should be positive";
//         EXPECT_GT(box[2], 0) << "Box length c should be positive";
//         EXPECT_NEAR(box[3], 90.0, 1e-6) << "Box angle alpha should be 90 degrees";
//         EXPECT_NEAR(box[4], 90.0, 1e-6) << "Box angle beta should be 90 degrees";
//         EXPECT_NEAR(box[5], 90.0, 1e-6) << "Box angle gamma should be 90 degrees";
//     }

//     // Compare missing topology information
//     auto psf_missing = psf_parser.get_missing_topology_info(pdb_atoms);
//     auto itp_missing = itp_parser.get_missing_topology_info(pdb_atoms);
    
//     EXPECT_EQ(psf_missing.size(), itp_missing.size()) 
//         << "Number of residues with missing topology mismatch";
    
//     for (const auto& [res, atoms] : psf_missing) {
//         EXPECT_TRUE(itp_missing.find(res) != itp_missing.end())
//             << "Residue " << res << " missing in ITP but not in PSF";
//         if (itp_missing.find(res) != itp_missing.end()) {
//             EXPECT_EQ(atoms, itp_missing[res])
//                 << "Missing atoms mismatch for residue " << res;
//         }
//     }
// }
