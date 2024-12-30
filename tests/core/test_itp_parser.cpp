// tests/core/test_itp_parser.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/io/itp_parser.hpp"
#include "config.hpp"
#include <cmath>

using namespace pygcmc::core::io;

TEST(ITPParserTest, ParseBENXFile) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/benx.itp";
    ASSERT_TRUE(parser.parse(itp_file)) << "Failed to parse BENX ITP file";

    // Test getting atom properties for BENX residue
    double charge, mass;
    
    // Test CG atom
    ASSERT_TRUE(parser.get_atom_properties("BENX", "CG", charge, mass)) 
        << "Failed to get properties for BENX CG atom";
    EXPECT_DOUBLE_EQ(charge, -0.115) << "Incorrect charge for BENX CG atom";
    EXPECT_DOUBLE_EQ(mass, 12.011) << "Incorrect mass for BENX CG atom";

    // Test HG atom
    ASSERT_TRUE(parser.get_atom_properties("BENX", "HG", charge, mass))
        << "Failed to get properties for BENX HG atom";
    EXPECT_DOUBLE_EQ(charge, 0.115) << "Incorrect charge for BENX HG atom";
    EXPECT_DOUBLE_EQ(mass, 1.008) << "Incorrect mass for BENX HG atom";

    // Test CD1 atom
    ASSERT_TRUE(parser.get_atom_properties("BENX", "CD1", charge, mass))
        << "Failed to get properties for BENX CD1 atom";
    EXPECT_DOUBLE_EQ(charge, -0.115) << "Incorrect charge for BENX CD1 atom";
    EXPECT_DOUBLE_EQ(mass, 12.011) << "Incorrect mass for BENX CD1 atom";

    // Test non-existent atom
    ASSERT_FALSE(parser.get_atom_properties("BENX", "XXX", charge, mass))
        << "Should not find properties for non-existent atom XXX";
}

TEST(ITPParserTest, UpdateBENXAtoms) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/benx.itp";
    ASSERT_TRUE(parser.parse(itp_file)) << "Failed to parse BENX ITP file";

    // Create test PDB atoms for BENX
    std::vector<PDBAtom> pdb_atoms;
    
    // Add BENX atoms with detailed comments
    PDBAtom cg(1, "CG", "BENX", 1);   // Carbon atom in benzene ring
    PDBAtom hg(2, "HG", "BENX", 1);   // Hydrogen attached to CG
    PDBAtom cd1(3, "CD1", "BENX", 1); // Carbon atom adjacent to CG
    
    pdb_atoms.push_back(cg);
    pdb_atoms.push_back(hg);
    pdb_atoms.push_back(cd1);

    // Update atoms and verify count
    int updated = parser.update_pdb_atoms(pdb_atoms);
    EXPECT_EQ(updated, 3) << "Expected 3 atoms to be updated";

    // Verify updated atoms with detailed messages
    // CG atom verification
    EXPECT_EQ(pdb_atoms[0].topo_type, "CG2R61") 
        << "Incorrect topology type for CG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_charge, -0.115) 
        << "Incorrect charge for CG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_mass, 12.011) 
        << "Incorrect mass for CG atom";

    // HG atom verification
    EXPECT_EQ(pdb_atoms[1].topo_type, "HGR61") 
        << "Incorrect topology type for HG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_charge, 0.115) 
        << "Incorrect charge for HG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_mass, 1.008) 
        << "Incorrect mass for HG atom";

    // CD1 atom verification
    EXPECT_EQ(pdb_atoms[2].topo_type, "CG2R61") 
        << "Incorrect topology type for CD1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[2].topo_charge, -0.115) 
        << "Incorrect charge for CD1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[2].topo_mass, 12.011) 
        << "Incorrect mass for CD1 atom";
}

TEST(ITPParserTest, ParseSOLFile) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/sol.itp";
    ASSERT_TRUE(parser.parse(itp_file)) << "Failed to parse SOL ITP file";

    // Test getting atom properties for SOL residue
    double charge, mass;
    
    // Test OW atom (Oxygen in water)
    ASSERT_TRUE(parser.get_atom_properties("SOL", "OW", charge, mass))
        << "Failed to get properties for SOL OW atom";
    EXPECT_DOUBLE_EQ(charge, -0.834) << "Incorrect charge for SOL OW atom";
    EXPECT_DOUBLE_EQ(mass, 15.9994) << "Incorrect mass for SOL OW atom";

    // Test HW1 atom (First Hydrogen in water)
    ASSERT_TRUE(parser.get_atom_properties("SOL", "HW1", charge, mass))
        << "Failed to get properties for SOL HW1 atom";
    EXPECT_DOUBLE_EQ(charge, 0.417) << "Incorrect charge for SOL HW1 atom";
    EXPECT_DOUBLE_EQ(mass, 1.008) << "Incorrect mass for SOL HW1 atom";

    // Test HW2 atom (Second Hydrogen in water)
    ASSERT_TRUE(parser.get_atom_properties("SOL", "HW2", charge, mass))
        << "Failed to get properties for SOL HW2 atom";
    EXPECT_DOUBLE_EQ(charge, 0.417) << "Incorrect charge for SOL HW2 atom";
    EXPECT_DOUBLE_EQ(mass, 1.008) << "Incorrect mass for SOL HW2 atom";

    // Verify charge neutrality of water molecule
    double total_charge = -0.834 + 0.417 + 0.417;
    EXPECT_NEAR(total_charge, 0.0, 1e-6) 
        << "Water molecule should have neutral total charge";
}

TEST(ITPParserTest, ParsePRPFile) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/prpx.itp";
    ASSERT_TRUE(parser.parse(itp_file)) << "Failed to parse PRPX ITP file";

    // Test getting atom properties for PRPX residue
    double charge, mass;
    
    // Test C1 atom (First Carbon in propane)
    ASSERT_TRUE(parser.get_atom_properties("PRPX", "C1", charge, mass))
        << "Failed to get properties for PRPX C1 atom";
    EXPECT_DOUBLE_EQ(charge, -0.27) << "Incorrect charge for PRPX C1 atom";
    EXPECT_DOUBLE_EQ(mass, 12.011) << "Incorrect mass for PRPX C1 atom";

    // Test H11 atom (First Hydrogen on C1)
    ASSERT_TRUE(parser.get_atom_properties("PRPX", "H11", charge, mass))
        << "Failed to get properties for PRPX H11 atom";
    EXPECT_DOUBLE_EQ(charge, 0.09) << "Incorrect charge for PRPX H11 atom";
    EXPECT_DOUBLE_EQ(mass, 1.008) << "Incorrect mass for PRPX H11 atom";

    // Test C2 atom (Second Carbon in propane)
    ASSERT_TRUE(parser.get_atom_properties("PRPX", "C2", charge, mass))
        << "Failed to get properties for PRPX C2 atom";
    EXPECT_DOUBLE_EQ(charge, -0.18) << "Incorrect charge for PRPX C2 atom";
    EXPECT_DOUBLE_EQ(mass, 12.011) << "Incorrect mass for PRPX C2 atom";
}

TEST(ITPParserTest, MissingTopologyInfo) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/benx.itp";
    ASSERT_TRUE(parser.parse(itp_file)) << "Failed to parse BENX ITP file";

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
    
    // Verify missing topology information
    EXPECT_TRUE(missing_info["BENX"].find("XXX") != missing_info["BENX"].end())
        << "XXX atom should be reported as missing for BENX";
    EXPECT_TRUE(missing_info["SOL"].find("OW") != missing_info["SOL"].end())
        << "OW atom should be reported as missing for SOL";
    EXPECT_TRUE(missing_info["BENX"].find("CG") == missing_info["BENX"].end())
        << "CG atom should not be reported as missing for BENX";
}

TEST(ITPParserTest, MultipleResidues) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/benx.itp";
    ASSERT_TRUE(parser.parse(itp_file)) << "Failed to parse BENX ITP file";

    // Create test PDB atoms for multiple BENX residues
    std::vector<PDBAtom> pdb_atoms;
    
    // First BENX residue (sequence number 1)
    PDBAtom cg1(1, "CG", "BENX", 1);
    PDBAtom hg1(2, "HG", "BENX", 1);
    
    // Second BENX residue (sequence number 2)
    PDBAtom cg2(3, "CG", "BENX", 2);
    PDBAtom hg2(4, "HG", "BENX", 2);
    
    pdb_atoms.push_back(cg1);
    pdb_atoms.push_back(hg1);
    pdb_atoms.push_back(cg2);
    pdb_atoms.push_back(hg2);

    // Update atoms and verify count
    int updated = parser.update_pdb_atoms(pdb_atoms);
    EXPECT_EQ(updated, 4) << "Expected 4 atoms to be updated";

    // Verify updated atoms for first residue
    EXPECT_EQ(pdb_atoms[0].topo_type, "CG2R61") 
        << "Incorrect topology type for first residue CG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_charge, -0.115) 
        << "Incorrect charge for first residue CG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_mass, 12.011) 
        << "Incorrect mass for first residue CG atom";

    EXPECT_EQ(pdb_atoms[1].topo_type, "HGR61") 
        << "Incorrect topology type for first residue HG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_charge, 0.115) 
        << "Incorrect charge for first residue HG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_mass, 1.008) 
        << "Incorrect mass for first residue HG atom";

    // Verify updated atoms for second residue
    EXPECT_EQ(pdb_atoms[2].topo_type, "CG2R61") 
        << "Incorrect topology type for second residue CG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[2].topo_charge, -0.115) 
        << "Incorrect charge for second residue CG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[2].topo_mass, 12.011) 
        << "Incorrect mass for second residue CG atom";

    EXPECT_EQ(pdb_atoms[3].topo_type, "HGR61") 
        << "Incorrect topology type for second residue HG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[3].topo_charge, 0.115) 
        << "Incorrect charge for second residue HG atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[3].topo_mass, 1.008) 
        << "Incorrect mass for second residue HG atom";
}

TEST(ITPParserTest, MultiplePRPXResidues) {
    ITPParser parser;
    std::string itp_file = std::string(ITP_DATA_DIR) + "/prpx.itp";
    ASSERT_TRUE(parser.parse(itp_file)) << "Failed to parse PRPX ITP file";

    // Create test PDB atoms for multiple PRPX residues
    std::vector<PDBAtom> pdb_atoms;
    
    // First PRPX residue (sequence number 1)
    PDBAtom c1_1(1, "C1", "PRPX", 1);
    PDBAtom h11_1(2, "H11", "PRPX", 1);
    PDBAtom h12_1(3, "H12", "PRPX", 1);
    
    // Second PRPX residue (sequence number 2)
    PDBAtom c1_2(4, "C1", "PRPX", 2);
    PDBAtom h11_2(5, "H11", "PRPX", 2);
    PDBAtom h12_2(6, "H12", "PRPX", 2);
    
    // Third PRPX residue (sequence number 3)
    PDBAtom c1_3(7, "C1", "PRPX", 3);
    PDBAtom h11_3(8, "H11", "PRPX", 3);
    PDBAtom h12_3(9, "H12", "PRPX", 3);
    
    pdb_atoms.push_back(c1_1);
    pdb_atoms.push_back(h11_1);
    pdb_atoms.push_back(h12_1);
    pdb_atoms.push_back(c1_2);
    pdb_atoms.push_back(h11_2);
    pdb_atoms.push_back(h12_2);
    pdb_atoms.push_back(c1_3);
    pdb_atoms.push_back(h11_3);
    pdb_atoms.push_back(h12_3);

    // Update atoms and verify count
    int updated = parser.update_pdb_atoms(pdb_atoms);
    EXPECT_EQ(updated, 9) << "Expected 9 atoms to be updated";

    // Verify updated atoms for first residue
    EXPECT_EQ(pdb_atoms[0].topo_type, "CG331") 
        << "Incorrect topology type for first residue C1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_charge, -0.27) 
        << "Incorrect charge for first residue C1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[0].topo_mass, 12.011) 
        << "Incorrect mass for first residue C1 atom";

    EXPECT_EQ(pdb_atoms[1].topo_type, "HGA3") 
        << "Incorrect topology type for first residue H11 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_charge, 0.09) 
        << "Incorrect charge for first residue H11 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[1].topo_mass, 1.008) 
        << "Incorrect mass for first residue H11 atom";

    // Verify updated atoms for second residue
    EXPECT_EQ(pdb_atoms[3].topo_type, "CG331") 
        << "Incorrect topology type for second residue C1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[3].topo_charge, -0.27) 
        << "Incorrect charge for second residue C1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[3].topo_mass, 12.011) 
        << "Incorrect mass for second residue C1 atom";

    // Verify updated atoms for third residue
    EXPECT_EQ(pdb_atoms[6].topo_type, "CG331") 
        << "Incorrect topology type for third residue C1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[6].topo_charge, -0.27) 
        << "Incorrect charge for third residue C1 atom";
    EXPECT_DOUBLE_EQ(pdb_atoms[6].topo_mass, 12.011) 
        << "Incorrect mass for third residue C1 atom";

    // Verify charge conservation for each residue
    double total_charge_res1 = pdb_atoms[0].topo_charge + 
                              pdb_atoms[1].topo_charge + 
                              pdb_atoms[2].topo_charge;
    EXPECT_NEAR(total_charge_res1, -0.09, 1e-6) 
        << "First PRPX residue should have expected total charge";

    double total_charge_res2 = pdb_atoms[3].topo_charge + 
                              pdb_atoms[4].topo_charge + 
                              pdb_atoms[5].topo_charge;
    EXPECT_NEAR(total_charge_res2, -0.09, 1e-6) 
        << "Second PRPX residue should have expected total charge";

    double total_charge_res3 = pdb_atoms[6].topo_charge + 
                              pdb_atoms[7].topo_charge + 
                              pdb_atoms[8].topo_charge;
    EXPECT_NEAR(total_charge_res3, -0.09, 1e-6) 
        << "Third PRPX residue should have expected total charge";
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}