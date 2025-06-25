# src/bindings/CMakeLists.txt

# First try to find system-installed pybind11
find_package(pybind11 CONFIG)

if(NOT pybind11_FOUND)
    message(STATUS "System pybind11 not found, downloading from GitHub...")
    # Download pybind11 if not found in system
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        # GIT_TAG        v2.11.1  # specify a version/tag
    )
    FetchContent_MakeAvailable(pybind11)
else()
    message(STATUS "Using system-installed pybind11")
endif()

# Automatically discover all source files
file(GLOB_RECURSE BINDING_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

# Automatically discover all header files (for IDE display)
file(GLOB_RECURSE BINDING_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx"
)

# Create Python module
pybind11_add_module(pygcmc 
    ${BINDING_SOURCES}
    ${BINDING_HEADERS}
)

# Link against other libraries
target_link_libraries(pygcmc
    PRIVATE
        pygcmc_model
        pygcmc_io
        pygcmc_system
        pygcmc_simulation
)

# Set position independent code (following system/io pattern)
set_property(TARGET pygcmc PROPERTY POSITION_INDEPENDENT_CODE ON)

# Set output directory to build root
set_target_properties(pygcmc PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
)

# Set include directories (following system/io pattern)
target_include_directories(pygcmc
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
) 