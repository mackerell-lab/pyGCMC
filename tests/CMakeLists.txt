# tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# 使用 FetchContent 获取 Catch2 v2.13.10
include(FetchContent)
FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.10
)
FetchContent_MakeAvailable(catch2)

# 添加测试可执行文件
add_executable(test_system core/test_system.cpp)

# 链接 Catch2 和核心库
target_link_libraries(test_system PRIVATE Catch2::Catch2 pygcmc_core)

# 启用 CTest
include(CTest)

# 注册测试
add_test(NAME test_system COMMAND test_system)
