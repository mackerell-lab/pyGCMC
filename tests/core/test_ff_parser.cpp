// tests/core/test_ff_parser.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/io/ff_parser.hpp"
#include <string>
#include <cmath>

namespace pygcmc {
namespace core {
namespace io {
namespace test {

// Global variable for parameter file path
std::string g_prm_file;

// Parse command line arguments
class Environment : public ::testing::Environment {
public:
    ~Environment() override {}

    void SetUp() override {
        // Get command line arguments
        const auto& args = ::testing::internal::GetArgvs();
        
        // Look for --prm_file argument
        for (const auto& arg : args) {
            if (arg.rfind("--prm_file=", 0) == 0) {
                g_prm_file = arg.substr(11);
                break;
            }
        }

        if (g_prm_file.empty()) {
            std::cerr << "Error: --prm_file argument is required" << std::endl;
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
    ASSERT_TRUE(parser_->parse(g_prm_file));
    const auto& params = parser_->get_nonbonded_params();

    // Test some atom types from par_all36_cgenff.prm
    // CG2R61 (aromatic carbon)
    {
        auto it = params.find("CG2R61");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.0700, 1e-4);
        EXPECT_NEAR(it->second.rmin, 1.9924, 1e-4);
    }

    // NG2R61 (aromatic nitrogen)
    {
        auto it = params.find("NG2R61");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.2000, 1e-4);
        EXPECT_NEAR(it->second.rmin, 1.8500, 1e-4);
    }

    // OG2D1 (carbonyl oxygen)
    {
        auto it = params.find("OG2D1");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.1200, 1e-4);
        EXPECT_NEAR(it->second.rmin, 1.7000, 1e-4);
    }

    // HGA1 (nonpolar H)
    {
        auto it = params.find("HGA1");
        ASSERT_NE(it, params.end());
        EXPECT_NEAR(it->second.epsilon, 0.0450, 1e-4);
        EXPECT_NEAR(it->second.rmin, 1.3400, 1e-4);
    }
}

// Test NBFIX parameters
TEST_F(TestFFParser, NBFIXParameters) {
    ASSERT_TRUE(parser_->parse(g_prm_file));
    const auto& nbfix = parser_->get_nbfix_params();

    // Test NBFIX pairs from par_all36_cgenff.prm
    // NG2S2-CLGR1 pair
    {
        auto key = std::make_pair("NG2S2", "CLGR1");
        auto it = nbfix.find(key);
        ASSERT_NE(it, nbfix.end());
        EXPECT_NEAR(it->second.epsilon, 0.40, 1e-4);
        EXPECT_NEAR(it->second.rmin, 3.88, 1e-4);

        // Test symmetry (CLGR1-NG2S2 should have same parameters)
        auto key_rev = std::make_pair("CLGR1", "NG2S2");
        auto it_rev = nbfix.find(key_rev);
        ASSERT_NE(it_rev, nbfix.end());
        EXPECT_NEAR(it_rev->second.epsilon, 0.40, 1e-4);
        EXPECT_NEAR(it_rev->second.rmin, 3.88, 1e-4);
    }

    // NG2S1-CLGR1 pair
    {
        auto key = std::make_pair("NG2S1", "CLGR1");
        auto it = nbfix.find(key);
        ASSERT_NE(it, nbfix.end());
        EXPECT_NEAR(it->second.epsilon, 0.40, 1e-4);
        EXPECT_NEAR(it->second.rmin, 3.88, 1e-4);
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

} // namespace test
} // namespace io
} // namespace core
} // namespace pygcmc

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new pygcmc::core::io::test::Environment);
    return RUN_ALL_TESTS();
}