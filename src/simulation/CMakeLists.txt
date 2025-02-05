# src/simulation/CMakeLists.txt

# Create simulation library
add_library(pygcmc_simulation STATIC
    simulation.hpp
    simulation.cpp
    # Future source files can be added here
)

# Set position independent code
set_property(TARGET pygcmc_simulation PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add include directories for this target
target_include_directories(pygcmc_simulation
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
)

# Link dependencies (simulation now depends on platform and model libraries)
target_link_libraries(pygcmc_simulation
    PUBLIC
        pygcmc_platform
        pygcmc_model
)

# Set properties
set_target_properties(pygcmc_simulation PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
) 