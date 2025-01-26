cmake_minimum_required(VERSION 3.15)
project(pygcmc_dev VERSION 1.0)

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
add_subdirectory(src/data)
add_subdirectory(src/forcefield)

# 2. Libraries depending on data
add_subdirectory(src/io)

# 3. System (depends on data and forcefield)
add_subdirectory(src/system)

# 4. Simulation (depends on system)
add_subdirectory(src/simulation)

# 5. Python bindings (depends on all other modules)
add_subdirectory(src/bindings)

# Add tests
add_subdirectory(tests)
