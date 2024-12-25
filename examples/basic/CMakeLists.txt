# examples/basic/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# 添加可执行文件
add_executable(run_system run_system.cpp)

# 链接核心库
target_link_libraries(run_system PRIVATE pygcmc_core)

# 包含目录
target_include_directories(run_system PRIVATE
    ../../modules/core/include
)
