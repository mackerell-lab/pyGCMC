# Create system library
add_library(pygcmc_system
    system.cpp
)

# Link dependencies - system depends on model
target_link_libraries(pygcmc_system
    PUBLIC 
        pygcmc_model      # system's interface exposes model types
)

# Set include directories for this target
target_include_directories(pygcmc_system
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

# Set properties
set_target_properties(pygcmc_system PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
) 