# modules/bindings/CMakeLists.txt

pybind11_add_module(pygcmc src/bindings.cpp)

target_link_libraries(pygcmc PRIVATE pygcmc_core)

target_include_directories(pygcmc
    PRIVATE
        ${CMAKE_SOURCE_DIR}/modules/core/include
)

# Set output directory
set_target_properties(pygcmc PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/modules/bindings
)
