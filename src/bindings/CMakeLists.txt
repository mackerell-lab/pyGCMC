# ./src/bindings/CMakeLists.txt

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

# Create Python module
pybind11_add_module(pygcmc 
    bindings.cpp
    model_bindings.cpp
    io_bindings.cpp
)

# Link against other libraries
target_link_libraries(pygcmc
    PRIVATE
        pygcmc_model
        pygcmc_io
)

# Set output directory to build root
set_target_properties(pygcmc PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Set include directories
target_include_directories(pygcmc
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
) 