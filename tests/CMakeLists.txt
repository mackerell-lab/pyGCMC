# tests/CMakeLists.txt

# 添加 PDBParser 测试
add_executable(test_pdb_parser core/test_pdb_parser.cpp)

# 链接 pygcmc_core 和 GTest
target_link_libraries(test_pdb_parser PRIVATE pygcmc_core GTest::gtest GTest::gtest_main)

# 添加 System 测试
add_executable(test_system core/test_system.cpp)

# 链接 pygcmc_core 和 GTest
target_link_libraries(test_system PRIVATE pygcmc_core GTest::gtest GTest::gtest_main)

# 注册测试
add_test(NAME TestPDBParser COMMAND test_pdb_parser)
add_test(NAME TestSystem COMMAND test_system)
