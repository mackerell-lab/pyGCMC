# src/model/CMakeLists.txt

# Automatically discover all header files (for IDE display and dependency tracking)
file(GLOB_RECURSE MODEL_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx"
)

# Create header-only model library
add_library(pygcmc_model INTERFACE)

# Add header files to target (for IDE organization)
target_sources(pygcmc_model INTERFACE ${MODEL_HEADERS})

# Add include directories
target_include_directories(pygcmc_model
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set properties
set_target_properties(pygcmc_model PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
) 