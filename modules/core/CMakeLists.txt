# modules/core/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# 定义核心库
add_library(pygcmc_core STATIC
    src/system.cpp
)

# 包含头文件目录
target_include_directories(pygcmc_core PUBLIC
    include
)

# 如果需要链接其他库，可以在这里添加
# target_link_libraries(pygcmc_core PUBLIC other_library)
