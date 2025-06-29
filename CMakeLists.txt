# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(pygcmc VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not set
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Set global include directory
include_directories(${PROJECT_SOURCE_DIR}/src)

# ==============================================================================
# UTILITY FUNCTIONS TO AVOID REPETITION
# ==============================================================================

# Function to create a standard library with common settings
function(add_pygcmc_library target_name source_dir library_type)
    # Auto-discover source files
    file(GLOB_RECURSE LIB_SOURCES
        "${source_dir}/*.cpp"
        "${source_dir}/*.cc"
        "${source_dir}/*.cxx"
    )
    
    # Auto-discover header files
    file(GLOB_RECURSE LIB_HEADERS
        "${source_dir}/*.hpp"
        "${source_dir}/*.h"
        "${source_dir}/*.hxx"
    )
    
    # Create library
    if(library_type STREQUAL "INTERFACE")
        add_library(${target_name} INTERFACE)
        target_sources(${target_name} INTERFACE ${LIB_HEADERS})
        target_include_directories(${target_name}
            INTERFACE
                ${source_dir}/..
                ${CMAKE_SOURCE_DIR}/src
                ${source_dir}
        )
    else()
        add_library(${target_name} ${library_type}
            ${LIB_SOURCES}
            ${LIB_HEADERS}
        )
        
        # Set position independent code
        set_property(TARGET ${target_name} PROPERTY POSITION_INDEPENDENT_CODE ON)
        
        # Add include directories
        target_include_directories(${target_name}
            PUBLIC
                ${source_dir}/..
                ${CMAKE_SOURCE_DIR}/src
        )
    endif()
    
    # Set common properties
    set_target_properties(${target_name} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        LINKER_LANGUAGE CXX
    )
endfunction()

# ==============================================================================
# BUILD LIBRARIES IN DEPENDENCY ORDER
# ==============================================================================

# 1. Model library (header-only, no dependencies)
add_pygcmc_library(pygcmc_model "${CMAKE_CURRENT_SOURCE_DIR}/src/model" INTERFACE)

# 2. IO library (depends on model)
add_pygcmc_library(pygcmc_io "${CMAKE_CURRENT_SOURCE_DIR}/src/io" STATIC)
target_link_libraries(pygcmc_io PUBLIC pygcmc_model)

# 3. System library (depends on model)
add_pygcmc_library(pygcmc_system "${CMAKE_CURRENT_SOURCE_DIR}/src/system" STATIC)
target_link_libraries(pygcmc_system PUBLIC pygcmc_model)

# 4. Platform library (depends on model)
add_pygcmc_library(pygcmc_platform "${CMAKE_CURRENT_SOURCE_DIR}/src/platform" STATIC)
target_link_libraries(pygcmc_platform PUBLIC pygcmc_model)

# 5. Simulation library (depends on platform and model)
add_library(pygcmc_simulation STATIC
    src/simulation/simulation.hpp
    src/simulation/simulation.cpp
)
set_property(TARGET pygcmc_simulation PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(pygcmc_simulation
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(pygcmc_simulation
    PUBLIC
        pygcmc_platform
        pygcmc_model
)
set_target_properties(pygcmc_simulation PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
)

# ==============================================================================
# PYTHON BINDINGS
# ==============================================================================

# Find or download pybind11
find_package(pybind11 CONFIG)
if(NOT pybind11_FOUND)
    message(STATUS "System pybind11 not found, downloading from GitHub...")
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
    )
    FetchContent_MakeAvailable(pybind11)
else()
    message(STATUS "Using system-installed pybind11")
endif()

# Auto-discover binding source files
file(GLOB_RECURSE BINDING_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.cxx"
)

# Auto-discover binding header files
file(GLOB_RECURSE BINDING_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.hxx"
)

# Create Python module
pybind11_add_module(pygcmc 
    ${BINDING_SOURCES}
    ${BINDING_HEADERS}
)

# Link against all libraries
target_link_libraries(pygcmc
    PRIVATE
        pygcmc_model
        pygcmc_io
        pygcmc_system
        pygcmc_platform
        pygcmc_simulation
)

# Set binding properties
set_property(TARGET pygcmc PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(pygcmc PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
)

target_include_directories(pygcmc
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

# ==============================================================================
# TESTING
# ==============================================================================

# Optional: Enable testing
option(BUILD_TESTING "Build tests" OFF)
if(BUILD_TESTING)
    enable_testing()
    # Uncomment when tests are ready
    # add_subdirectory(tests)
endif()