# Add Python bindings
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.11.1  # specify a version/tag
)
FetchContent_MakeAvailable(pybind11)

# Create Python module
pybind11_add_module(pygcmc_python 
    io_bindings.cpp
)

# Link against other libraries
target_link_libraries(pygcmc_python
    PRIVATE
        pygcmc_model
        pygcmc_io
)

# Set include directories
target_include_directories(pygcmc_python
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${pybind11_INCLUDE_DIRS}  # Add pybind11 include directories
) 