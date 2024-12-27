// tests/core/test_pdb_parser.cpp

#include "gtest/gtest.h"
#include "pygcmc/core/io/pdb_parser.hpp"

TEST(PDBParserTest, ParseValidFile) {
    std::string filename = "tests/data/test.pdb"; // 确保存在测试 PDB 文件
    auto parsed = pygcmc::core::io::PDBParser::parse(filename);
    
    EXPECT_FALSE(parsed.first.empty()) << "晶胞信息应非空";
    EXPECT_FALSE(parsed.second.empty()) << "原子列表应非空";
}

TEST(PDBParserTest, ParseInvalidFile) {
    std::string filename = "tests/data/nonexistent.pdb";
    EXPECT_THROW({
        pygcmc::core::io::PDBParser::parse(filename);
    }, std::runtime_error);
}
