# Create io library
add_library(pygcmc_io
    pdbParser.cpp
)

# Link dependencies - io depends on data
target_link_libraries(pygcmc_io
    PRIVATE pygcmc_data  # io implementation needs data, but users of io don't need data directly
)

# Set include directories for this target
target_include_directories(pygcmc_io
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

# Set properties
set_target_properties(pygcmc_io PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
) 