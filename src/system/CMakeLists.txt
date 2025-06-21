# src/system/CMakeLists.txt

# Automatically discover all source files
file(GLOB_RECURSE SYSTEM_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

# Automatically discover all header files (for IDE display)
file(GLOB_RECURSE SYSTEM_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx"
)

# Create system library
add_library(pygcmc_system STATIC
    ${SYSTEM_SOURCES}
    ${SYSTEM_HEADERS}
)

# Set position independent code
set_property(TARGET pygcmc_system PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add include directories
target_include_directories(pygcmc_system
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
)

# Link dependencies (system depends on the model library)
target_link_libraries(pygcmc_system
    PUBLIC
        pygcmc_model
)

# Set properties
set_target_properties(pygcmc_system PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
) 