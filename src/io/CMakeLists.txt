# src/io/CMakeLists.txt

# Automatically discover all source files
file(GLOB_RECURSE IO_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

# Automatically discover all header files (for IDE display)
file(GLOB_RECURSE IO_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx"
)

# Create IO library
add_library(pygcmc_io STATIC
    ${IO_SOURCES}
    ${IO_HEADERS}
)

# Set position independent code
set_property(TARGET pygcmc_io PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add include directories
target_include_directories(pygcmc_io
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
)

# Link dependencies
target_link_libraries(pygcmc_io
    PUBLIC
        pygcmc_model
)

# Set properties
set_target_properties(pygcmc_io PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
) 