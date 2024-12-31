// tests/core/test_ff_parser.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/io/ff_parser.hpp"
#include <string>
#include <cmath>

namespace pygcmc {
namespace core {
namespace io {
namespace test {

// Global variables for parameter file paths
std::string g_prm_file;
std::string g_prot_prm_file;

// Parse command line arguments
class Environment : public ::testing::Environment {
public:
    ~Environment() override {}

    void SetUp() override {
        // Get command line arguments
        const auto& args = ::testing::internal::GetArgvs();
        
        // Look for --prm_file and --prot_prm_file arguments
        for (const auto& arg : args) {
            if (arg.rfind("--prm_file=", 0) == 0) {
                g_prm_file = arg.substr(11);
            }
            else if (arg.rfind("--prot_prm_file=", 0) == 0) {
                g_prot_prm_file = arg.substr(16);
            }
        }

        if (g_prm_file.empty()) {
            std::cerr << "Error: --prm_file argument is required" << std::endl;
            exit(1);
        }
        if (g_prot_prm_file.empty()) {
            std::cerr << "Error: --prot_prm_file argument is required" << std::endl;
            exit(1);
        }
    }

    void TearDown() override {}
};

class TestFFParser : public ::testing::Test {
protected:
    void SetUp() override {
        parser_ = std::make_unique<FFParser>();
    }

    std::unique_ptr<FFParser> parser_;
};

// Test basic file parsing
TEST_F(TestFFParser, ParsePRMFile) {
    ASSERT_TRUE(parser_->parse(g_prm_file));
}

// Test nonbonded parameters
TEST_F(TestFFParser, NonbondedParameters) {
    // First parse the water parameters
    ASSERT_TRUE(parser_->parse(std::string(TEST_DATA_DIR) + "/toppar_water_ions.str"))
        << "Failed to parse water and ions STR file";
    
    const auto& params = parser_->get_nonbonded_params();

    // Test water oxygen parameters
    {
        auto it = params.find("OT");
        ASSERT_NE(it, params.end()) << "Failed to find OT atom type";
        EXPECT_NEAR(it->second.epsilon, 0.1521, 1e-4) << "Incorrect epsilon for OT atom";
        EXPECT_NEAR(it->second.rmin, 1.7682 * 2, 1e-4) << "Incorrect Rmin for OT atom";
    }

    // Test water hydrogen parameters
    {
        auto it = params.find("HT");
        ASSERT_NE(it, params.end()) << "Failed to find HT atom type";
        EXPECT_NEAR(it->second.epsilon, 0.046, 1e-4) << "Incorrect epsilon for HT atom";
        EXPECT_NEAR(it->second.rmin, 0.2245 * 2, 1e-4) << "Incorrect Rmin for HT atom";
    }
}

// Test NBFIX parameters
TEST_F(TestFFParser, NBFIXParameters) {
    // First parse the water parameters
    ASSERT_TRUE(parser_->parse(std::string(TEST_DATA_DIR) + "/toppar_water_ions.str"))
        << "Failed to parse water and ions STR file";
    
    const auto& nbfix = parser_->get_nbfix_params();

    // Test SOD-CLA pair
    {
        auto key = std::make_pair("SOD", "CLA");
        auto it = nbfix.find(key);
        ASSERT_NE(it, nbfix.end()) << "Failed to find SOD-CLA NBFIX parameters";
        EXPECT_NEAR(it->second.epsilon, 0.0839, 1e-4) << "Incorrect NBFIX epsilon for SOD-CLA";
        EXPECT_NEAR(it->second.rmin, 3.7310, 1e-4) << "Incorrect NBFIX Rmin for SOD-CLA";

        // Test symmetry (CLA-SOD should have same parameters)
        auto key_rev = std::make_pair("CLA", "SOD");
        auto it_rev = nbfix.find(key_rev);
        ASSERT_NE(it_rev, nbfix.end()) << "Failed to find CLA-SOD NBFIX parameters";
        EXPECT_NEAR(it_rev->second.epsilon, 0.0839, 1e-4) << "Incorrect NBFIX epsilon for CLA-SOD";
        EXPECT_NEAR(it_rev->second.rmin, 3.7310, 1e-4) << "Incorrect NBFIX Rmin for CLA-SOD";
    }

    // Test POT-CLA pair
    {
        auto key = std::make_pair("POT", "CLA");
        auto it = nbfix.find(key);
        ASSERT_NE(it, nbfix.end()) << "Failed to find POT-CLA NBFIX parameters";
        EXPECT_NEAR(it->second.epsilon, 0.1142, 1e-4) << "Incorrect NBFIX epsilon for POT-CLA";
        EXPECT_NEAR(it->second.rmin, 4.0810, 1e-4) << "Incorrect NBFIX Rmin for POT-CLA";
    }
}

// Test error handling
TEST_F(TestFFParser, ErrorHandling) {
    // Test non-existent file
    ASSERT_FALSE(parser_->parse("non_existent_file.prm"));

    // Test with valid file
    ASSERT_TRUE(parser_->parse(g_prm_file));
    const auto& params = parser_->get_nonbonded_params();

    // Test non-existent atom type
    EXPECT_EQ(params.find("NON_EXISTENT"), params.end());

    // Test non-existent NBFIX pair
    const auto& nbfix = parser_->get_nbfix_params();
    auto key = std::make_pair("NON_EXISTENT1", "NON_EXISTENT2");
    EXPECT_EQ(nbfix.find(key), nbfix.end());
}

// Test protein parameters
TEST_F(TestFFParser, ProteinParameters) {
    // First parse the water parameters
    ASSERT_TRUE(parser_->parse(std::string(TEST_DATA_DIR) + "/toppar_water_ions.str"))
        << "Failed to parse water parameters file";
    
    // Then parse the protein force field
    ASSERT_TRUE(parser_->parse(g_prot_prm_file))
        << "Failed to parse protein force field file";
    
    const auto& params = parser_->get_nonbonded_params();

    // Test NH3 (N-terminal nitrogen)
    {
        auto it = params.find("NH3");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.2000, 1e-4);
        EXPECT_NEAR(it->second.rmin, 3.7000, 1e-4);  // 1.85 * 2
    }

    // Test CT1 (alpha carbon)
    {
        auto it = params.find("CT1");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.0320, 1e-4);
        EXPECT_NEAR(it->second.rmin, 4.0000, 1e-4);  // 2.000 * 2
    }

    // Test O (carbonyl oxygen)
    {
        auto it = params.find("O");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.1200, 1e-4);
        EXPECT_NEAR(it->second.rmin, 3.4000, 1e-4);  // 1.700 * 2
    }

    // Test HA1 (alpha hydrogen)
    {
        auto it = params.find("HA1");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.0450, 1e-4);
        EXPECT_NEAR(it->second.rmin, 2.6800, 1e-4);  // 1.340 * 2
    }
}

// Test NBFIX parameters from protein force field
TEST_F(TestFFParser, ProteinNBFIXParameters) {
    // First parse the general force field
    ASSERT_TRUE(parser_->parse(g_prm_file))
        << "Failed to parse general force field file";
    
    // Then parse the protein force field
    ASSERT_TRUE(parser_->parse(g_prot_prm_file))
        << "Failed to parse protein force field file";
    
    const auto& nbfix = parser_->get_nbfix_params();

    // Test SOD-OC pair
    {
        auto key = std::make_pair("SOD", "OC");
        auto it = nbfix.find(key);
        ASSERT_NE(it, nbfix.end()) << "Failed to find SOD-OC NBFIX parameters";
        EXPECT_NEAR(it->second.epsilon, 0.07502, 1e-4) << "Incorrect NBFIX epsilon for SOD-OC";
        EXPECT_NEAR(it->second.rmin, 3.23, 1e-4) << "Incorrect NBFIX Rmin for SOD-OC";

        // Test symmetry (OC-SOD should have same parameters)
        auto key_rev = std::make_pair("OC", "SOD");
        auto it_rev = nbfix.find(key_rev);
        ASSERT_NE(it_rev, nbfix.end()) << "Failed to find OC-SOD NBFIX parameters";
        EXPECT_NEAR(it_rev->second.epsilon, 0.07502, 1e-4) << "Incorrect NBFIX epsilon for OC-SOD";
        EXPECT_NEAR(it_rev->second.rmin, 3.23, 1e-4) << "Incorrect NBFIX Rmin for OC-SOD";
    }
}

// Test parsing of .str file format
TEST_F(TestFFParser, ParseSTRFile) {
    ASSERT_TRUE(parser_->parse(std::string(TEST_DATA_DIR) + "/toppar_water_ions.str"))
        << "Failed to parse water and ions STR file";
    
    const auto& params = parser_->get_nonbonded_params();
    const auto& nbfix = parser_->get_nbfix_params();

    // Test water oxygen parameters
    {
        auto it = params.find("OT");
        ASSERT_NE(it, params.end()) << "Failed to find OT atom type";
        EXPECT_NEAR(it->second.epsilon, 0.1521, 1e-4) << "Incorrect epsilon for OT atom";
        EXPECT_NEAR(it->second.rmin, 1.7682 * 2, 1e-4) << "Incorrect Rmin for OT atom";
    }

    // Test water hydrogen parameters
    {
        auto it = params.find("HT");
        ASSERT_NE(it, params.end()) << "Failed to find HT atom type";
        EXPECT_NEAR(it->second.epsilon, 0.046, 1e-4) << "Incorrect epsilon for HT atom";
        EXPECT_NEAR(it->second.rmin, 0.2245 * 2, 1e-4) << "Incorrect Rmin for HT atom";
    }

    // Test sodium parameters
    {
        auto it = params.find("SOD");
        ASSERT_NE(it, params.end()) << "Failed to find SOD atom type";
        EXPECT_NEAR(it->second.epsilon, 0.0469, 1e-4) << "Incorrect epsilon for SOD atom";
        EXPECT_NEAR(it->second.rmin, 1.41075 * 2, 1e-4) << "Incorrect Rmin for SOD atom";
    }

    // Test chloride parameters
    {
        auto it = params.find("CLA");
        ASSERT_NE(it, params.end()) << "Failed to find CLA atom type";
        EXPECT_NEAR(it->second.epsilon, 0.150, 1e-4) << "Incorrect epsilon for CLA atom";
        EXPECT_NEAR(it->second.rmin, 2.27 * 2, 1e-4) << "Incorrect Rmin for CLA atom";
    }
}

// Test water parameters
TEST_F(TestFFParser, WaterParameters) {
    // First parse the water parameters
    ASSERT_TRUE(parser_->parse(std::string(TEST_DATA_DIR) + "/toppar_water_ions.str"))
        << "Failed to parse water parameters file";
    
    const auto& params = parser_->get_nonbonded_params();

    // Test OT (TIP3P water oxygen)
    {
        auto it = params.find("OT");
        ASSERT_NE(it, params.end()) << "Failed to find OT parameters";
        EXPECT_NEAR(it->second.epsilon, 0.1521, 1e-4) << "Incorrect epsilon for OT";
        EXPECT_NEAR(it->second.rmin, 1.7682 * 2, 1e-4) << "Incorrect Rmin for OT";
    }

    // Test HT (TIP3P water hydrogen)
    {
        auto it = params.find("HT");
        ASSERT_NE(it, params.end()) << "Failed to find HT parameters";
        EXPECT_NEAR(it->second.epsilon, 0.046, 1e-4) << "Incorrect epsilon for HT";
        EXPECT_NEAR(it->second.rmin, 0.2245 * 2, 1e-4) << "Incorrect Rmin for HT";
    }
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