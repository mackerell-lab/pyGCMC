# tests/CMakeLists.txt

# 添加测试
enable_testing()

# 设置包含目录
include_directories(
    ${CMAKE_SOURCE_DIR}/modules/core/include
    ${CMAKE_BINARY_DIR}/tests
)

# 添加 TestPDBParser 测试
add_executable(test_pdb_parser core/test_pdb_parser.cpp)
target_link_libraries(test_pdb_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 添加 TestTopParser 测试
add_executable(test_top_parser core/test_top_parser.cpp)
target_link_libraries(test_top_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 设置 PDB 文件路径
set(TEST_PDB_FILE "${CMAKE_SOURCE_DIR}/tests/data/test.pdb")

# 添加测试并传递命令行参数
add_test(
    NAME TestPDBParser
    COMMAND test_pdb_parser --pdb_file=${TEST_PDB_FILE}
)

add_test(
    NAME TestTopParser
    COMMAND test_top_parser
)

# 添加 TestSystem 测试
add_executable(test_system core/test_system.cpp)
target_link_libraries(test_system PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_test(
    NAME TestSystem
    COMMAND test_system
)

# Configure config.hpp
configure_file(
    "${CMAKE_SOURCE_DIR}/tests/config.hpp.in"
    "${CMAKE_BINARY_DIR}/tests/config.hpp"
    @ONLY
)

# Add include directories for all test targets
target_include_directories(test_pdb_parser PRIVATE 
    ${CMAKE_BINARY_DIR}/tests
)

target_include_directories(test_top_parser PRIVATE 
    ${CMAKE_BINARY_DIR}/tests
)

target_include_directories(test_system PRIVATE 
    ${CMAKE_BINARY_DIR}/tests
)
