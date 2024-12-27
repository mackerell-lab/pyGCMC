# tests/CMakeLists.txt

# 添加测试
enable_testing()

# 添加 TestPDBParser 测试
add_executable(test_pdb_parser core/test_pdb_parser.cpp)
target_link_libraries(test_pdb_parser PRIVATE pygcmc_core gtest_main)

# 设置 PDB 文件路径
set(TEST_PDB_FILE "${CMAKE_SOURCE_DIR}/tests/data/test.pdb")

# 添加测试并传递命令行参数
add_test(
    NAME TestPDBParser
    COMMAND test_pdb_parser --pdb_file=${TEST_PDB_FILE}
)

# 添加其他测试
add_executable(test_system core/test_system.cpp)
target_link_libraries(test_system PRIVATE pygcmc_core gtest_main)

add_test(
    NAME TestSystem
    COMMAND test_system
)
