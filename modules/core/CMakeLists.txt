# modules/core/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# Define the core library
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
)

# Enable compiler warnings
if(MSVC)
    target_compile_options(pygcmc_core PRIVATE /W4)
else()
    target_compile_options(pygcmc_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Position Independent Code for potential shared libraries
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

# Link dependencies if any (e.g., OpenMP)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(pygcmc_core PUBLIC OpenMP::OpenMP_CXX)
endif()
