# Create IO library
add_library(pygcmc_io STATIC
    pdbParser.cpp
    psfParser.cpp
)

# Set position independent code
set_property(TARGET pygcmc_io PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add include directories
target_include_directories(pygcmc_io
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
)

# Link dependencies
target_link_libraries(pygcmc_io
    PUBLIC
        pygcmc_model
)

# Set properties
set_target_properties(pygcmc_io PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
) 