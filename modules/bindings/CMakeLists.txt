# modules/bindings/CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

# Define the bindings library
pybind11_add_module(pygcmc src/bindings.cpp)

# Link with core library and any other dependencies
target_link_libraries(pygcmc PRIVATE
    pygcmc_core
    pybind11::module
)

# Set include directories if needed
target_include_directories(pygcmc PRIVATE
    ${CMAKE_SOURCE_DIR}/modules/core/include
)

# Set output directory for Python modules
set_target_properties(pygcmc PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
)

# Installation rules
install(TARGETS pygcmc
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/python
)
