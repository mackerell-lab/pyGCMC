# Create header-only forcefield library
add_library(pygcmc_forcefield INTERFACE)

# Set include directories for this target
target_include_directories(pygcmc_forcefield
    INTERFACE
        ${CMAKE_SOURCE_DIR}/src
)

# Set properties
set_target_properties(pygcmc_forcefield PROPERTIES
    VERSION ${PROJECT_VERSION}
) 