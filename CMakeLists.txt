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

# Build in order of dependencies
# 1. Core libraries (no dependencies)
add_subdirectory(src/model)
# add_subdirectory(src/forcefield)  # Not implemented yet

# 2. Libraries depending on model
add_subdirectory(src/io)

# 3. System (depends on model and forcefield)
# add_subdirectory(src/system)      # Not implemented yet

# 4. Simulation (depends on system)
# add_subdirectory(src/simulation)  # Not implemented yet

# 5. Python bindings (depends on implemented modules)
add_subdirectory(src/bindings)

# Optional: Enable testing
option(BUILD_TESTING "Build tests" OFF)
if(BUILD_TESTING)
    enable_testing()
    # Uncomment when tests are ready
    # add_subdirectory(tests)
endif()
