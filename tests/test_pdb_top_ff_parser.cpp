TEST(PDBTOPFFParserTest, ReadPDBTOPFFTest) {
    std::string test_data_dir = "tests/data/";
    std::string pdb_file = test_data_dir + "test.pdb";
    std::string top_file = test_data_dir + "test.top";
    std::string cgenff_file = test_data_dir + "par_all36_cgenff.prm";
    std::string prot_file = test_data_dir + "par_all36m_prot.prm";

    PDBTOPFFParser parser;
    parser.read_pdb(pdb_file);
    parser.read_top(top_file);
    parser.read_ff(cgenff_file);
    parser.read_ff(prot_file);

    // Verify that both parameter files were read correctly
    ASSERT_TRUE(parser.has_ff_parameters());
    
    // Test some specific parameters from both files
    // CGenFF parameters
    auto cgenff_params = parser.get_ff_parameters();
    ASSERT_TRUE(cgenff_params.find("CG331") != cgenff_params.end());
    ASSERT_TRUE(cgenff_params.find("OG301") != cgenff_params.end());
    
    // Protein parameters
    ASSERT_TRUE(cgenff_params.find("NH1") != cgenff_params.end());
    ASSERT_TRUE(cgenff_params.find("CT1") != cgenff_params.end());
    ASSERT_TRUE(cgenff_params.find("CT2") != cgenff_params.end());
    
    // Test some specific bond parameters
    auto bonds = parser.get_bond_parameters();
    ASSERT_FALSE(bonds.empty());
    
    // Test some specific angle parameters
    auto angles = parser.get_angle_parameters();
    ASSERT_FALSE(angles.empty());
    
    // Test some specific dihedral parameters
    auto dihedrals = parser.get_dihedral_parameters();
    ASSERT_FALSE(dihedrals.empty());
    
    // Test some specific improper parameters
    auto impropers = parser.get_improper_parameters();
    ASSERT_FALSE(impropers.empty());
} 