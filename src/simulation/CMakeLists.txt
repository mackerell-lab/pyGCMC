# Create simulation library
add_library(pygcmc_simulation
    simulation.cpp
)

# Link dependencies - simulation depends on system
target_link_libraries(pygcmc_simulation
    PRIVATE pygcmc_system  # simulation implementation needs system
)

# Set include directories for this target
target_include_directories(pygcmc_simulation
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

# Set properties
set_target_properties(pygcmc_simulation PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
) 