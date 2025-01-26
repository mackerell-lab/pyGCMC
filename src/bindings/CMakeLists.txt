# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Create Python module
pybind11_add_module(pygcmc_python MODULE
    bindings.cpp
)

# Link all libraries that need Python bindings
target_link_libraries(pygcmc_python
    PRIVATE
        pygcmc_data
        pygcmc_forcefield
        pygcmc_io
        pygcmc_system
        pygcmc_simulation
)

# Set include directories
target_include_directories(pygcmc_python
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
) 