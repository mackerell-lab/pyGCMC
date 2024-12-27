# modules/core/CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Define source files groups
set(CORE_SOURCES
    src/system.cpp
    src/force_field.cpp
    src/box.cpp
)

# IO Parser source files
set(IO_SOURCES
    src/io/pdb_parser.cpp
    src/io/psf_parser.cpp
    src/io/top_parser.cpp
    src/io/itp_parser.cpp
    src/io/ff_parser.cpp
)

# GCMC source files
# 添加 GCMC 模块的源文件，如果有的话

# Optimization source files
# 添加 Optimization 模块的源文件，如果有的话

# Add core static library
add_library(pygcmc_core STATIC
    ${CORE_SOURCES}
    ${IO_SOURCES}
)

# Enable compiler warnings
if(MSVC)
    target_compile_options(pygcmc_core PRIVATE /W4)
else()
    target_compile_options(pygcmc_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Enable position independent code
set_target_properties(pygcmc_core PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    DEBUG_POSTFIX "d"
)

# Include directories
target_include_directories(pygcmc_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Create Python binding module
pybind11_add_module(pyGCMC MODULE
    src/bindings.cpp
)

# Link dependencies
target_link_libraries(pyGCMC 
    PRIVATE 
        pygcmc_core
        pybind11::module
)

# Set output directories
set_target_properties(pyGCMC PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
)

# Installation rules
install(TARGETS pygcmc_core pyGCMC
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)
