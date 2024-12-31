# tests/CMakeLists.txt

# 添加测试
enable_testing()

# 设置测试数据目录
set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
add_compile_definitions(TEST_DATA_DIR="${TEST_DATA_DIR}")

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

# 添加 TestPSFParser 测试
add_executable(test_psf_parser core/test_psf_parser.cpp)
target_link_libraries(test_psf_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 添加 TestITPParser 测试
add_executable(test_itp_parser core/test_itp_parser.cpp)
target_link_libraries(test_itp_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 添加 TestTopAndITP 测试
add_executable(test_top_and_itp core/test_top_and_itp_parser.cpp)
target_link_libraries(test_top_and_itp PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 添加 TestPSFAndITP 测试
add_executable(test_psf_and_itp core/test_psf_and_itp_parser.cpp)
target_link_libraries(test_psf_and_itp PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 添加 TestPDBTopFFParser 测试
add_executable(test_pdb_top_ff_parser core/test_pdb_top_ff_parser.cpp)
target_link_libraries(test_pdb_top_ff_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 添加 TestSystem 测试
add_executable(test_system core/test_system.cpp)
target_link_libraries(test_system PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# 添加所有测试到 CTest
add_test(NAME test_pdb_parser COMMAND test_pdb_parser)
add_test(NAME test_top_parser COMMAND test_top_parser)
add_test(NAME test_psf_parser COMMAND test_psf_parser)
add_test(NAME test_itp_parser COMMAND test_itp_parser)
add_test(NAME test_top_and_itp COMMAND test_top_and_itp)
add_test(NAME test_psf_and_itp COMMAND test_psf_and_itp)
add_test(NAME test_pdb_top_ff_parser COMMAND test_pdb_top_ff_parser)
add_test(NAME test_system COMMAND test_system)

# Configure config.hpp
configure_file(
    "${CMAKE_SOURCE_DIR}/tests/config.hpp.in"
    "${CMAKE_BINARY_DIR}/tests/config.hpp"
    @ONLY
)

# Add include directories for all test targets
target_include_directories(test_pdb_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_top_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_psf_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_itp_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_top_and_itp PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_psf_and_itp PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_pdb_top_ff_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_system PRIVATE ${CMAKE_BINARY_DIR}/tests)
