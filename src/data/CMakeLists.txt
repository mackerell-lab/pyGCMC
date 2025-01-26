# Create header-only data library
add_library(pygcmc_data INTERFACE)

# Set include directories for this target
target_include_directories(pygcmc_data
    INTERFACE
        ${CMAKE_SOURCE_DIR}/src
)

# Set properties
set_target_properties(pygcmc_data PROPERTIES
    VERSION ${PROJECT_VERSION}
) 