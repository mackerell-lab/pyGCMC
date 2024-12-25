# CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(pyGCMC LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加核心模块
add_subdirectory(modules/core)

# 添加测试
enable_testing()
add_subdirectory(tests)

# 添加示例代码
add_subdirectory(examples/basic)
