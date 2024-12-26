cmake_minimum_required(VERSION 3.15)
project(pyGCMC LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 尝试在系统中查找已安装的 pybind11
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND)
    message(STATUS "pybind11 未在系统中找到。将通过 FetchContent 下载 pybind11。")
    
    # 包含 FetchContent 模块
    include(FetchContent)
    
    # 声明并下载 pybind11
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG v2.10.1 # 您可以选择特定的版本标签
    )
    
    # 下载并添加 pybind11
    FetchContent_MakeAvailable(pybind11)
else()
    message(STATUS "已找到系统中的 pybind11: ${pybind11_DIR}")
endif()

# 添加核心模块
add_subdirectory(modules/core)

# 添加测试
enable_testing()
add_subdirectory(tests)

# 添加示例代码
add_subdirectory(examples/basic)
