// tests/core/test_pdb_parser.cpp

#include "gtest/gtest.h"
#include "pygcmc/core/io/pdb_parser.hpp"
#include <iostream>
#include <string>

// 使用一个全局变量来存储文件路径
std::string g_test_pdb_file;

// 主函数，用于初始化 Google Test 并解析命令行参数
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    
    // 解析命令行参数
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--pdb_file=") == 0) {
            g_test_pdb_file = arg.substr(11);
        }
    }
    
    return RUN_ALL_TESTS();
}

TEST(PDBParserTest, ParseValidFile) {
    ASSERT_FALSE(g_test_pdb_file.empty()) << "PDB 文件路径未提供";
    
    std::string filename = g_test_pdb_file;
    
    try {
        auto parsed = pygcmc::core::io::PDBParser::parse(filename);
        
        // 打印晶胞信息
        std::cout << "Crystal Parameters (a, b, c): ";
        for (const auto& param : parsed.first) {
            std::cout << param << " ";
        }
        std::cout << std::endl;
        
        // 打印原子信息
        std::cout << "Atoms:" << std::endl;
        for (const auto& atom : parsed.second) {
            std::cout << "Serial: " << atom.serial
                      << ", Name: " << atom.name
                      << ", Residue: " << atom.residue
                      << ", Sequence: " << atom.sequence
                      << ", Chain: " << atom.chain
                      << ", X: " << atom.x
                      << ", Y: " << atom.y
                      << ", Z: " << atom.z
                      << ", Occupancy: " << atom.occupancy
                      << ", Temp Factor: " << atom.temp_factor
                      << ", Element: " << atom.element
                      << ", Charge: " << atom.charge
                      << ", Type: " << atom.type
                      << std::endl;
        }
        
        // 断言解析结果不为空
        EXPECT_FALSE(parsed.first.empty()) << "晶胞信息应非空";
        EXPECT_FALSE(parsed.second.empty()) << "原子列表应非空";
    } catch (const std::exception& e) {
        FAIL() << "解析 PDB 文件时抛出异常: " << e.what();
    }
}

TEST(PDBParserTest, ParseInvalidFile) {
    std::string filename = "../tests/data/nonexistent.pdb"; // 仍使用相对路径，确保路径无效
    EXPECT_THROW({
        pygcmc::core::io::PDBParser::parse(filename);
    }, pygcmc::core::io::ParserError);
}
