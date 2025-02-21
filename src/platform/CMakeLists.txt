# src/platform/CMakeLists.txt

# Create platform library
add_library(pygcmc_platform STATIC
    platform.hpp
    cpu/energy.hpp
    cpu/energyCommon.hpp
    cpu/energyCommon.cpp
    cpu/energyDirect.hpp
    cpu/energyEwald.hpp
    cpu/energyEwald.cpp
    cpu/energyDirect.cpp
    # Future source files like cpuPlatform.cpp, cudaPlatform.cu, etc. can be added here
)

# Set position independent code
set_property(TARGET pygcmc_platform PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add include directories
target_include_directories(pygcmc_platform
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/src
)

# Link dependencies (for example, platform depends on the model library)
target_link_libraries(pygcmc_platform
    PUBLIC
        pygcmc_model
)

# Set properties
set_target_properties(pygcmc_platform PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    LINKER_LANGUAGE CXX
)
