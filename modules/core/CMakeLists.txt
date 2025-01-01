# modules/core/CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

# 定义核心库
add_library(pygcmc_core STATIC
    src/system.cpp
    src/force_field.cpp
    src/box.cpp
    src/io/pdb_parser.cpp
    src/io/psf_parser.cpp
    src/io/top_parser.cpp
    src/io/itp_parser.cpp
    src/io/ff_parser.cpp
    src/particles.cpp
    src/project.cpp
    src/structure.cpp
    src/project_atom.cpp
    src/project_residue.cpp
)

# 启用编译器警告
if(MSVC)
    target_compile_options(pygcmc_core PRIVATE /W4)
else()
    target_compile_options(pygcmc_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

# 设置位置独立代码（用于共享库）
set_target_properties(pygcmc_core PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    DEBUG_POSTFIX "d"
)

# 包含目录
target_include_directories(pygcmc_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# 链接依赖（例如 OpenMP）
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(pygcmc_core PUBLIC OpenMP::OpenMP_CXX)
endif()
