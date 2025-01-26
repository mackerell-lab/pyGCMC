# Add Python bindings
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.11.1  # specify a version/tag
)
FetchContent_MakeAvailable(pybind11)

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