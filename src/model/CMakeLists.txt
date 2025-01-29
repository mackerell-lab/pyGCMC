# src/model/CMakeLists.txt

# Create header-only model library
add_library(pygcmc_model INTERFACE)

# Add include directories
target_include_directories(pygcmc_model
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
)

# Set properties
set_target_properties(pygcmc_model PROPERTIES
    LINKER_LANGUAGE CXX
) 