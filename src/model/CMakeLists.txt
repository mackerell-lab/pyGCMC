# Create header-only model library
add_library(pygcmc_model INTERFACE)

# Add include directories
target_include_directories(pygcmc_model
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(pygcmc_model PROPERTIES
    LINKER_LANGUAGE CXX
) 