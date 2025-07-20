cmake_minimum_required(VERSION 3.10)

# Create benchmark executable
add_executable(benchmark_drude_vs_standard benchmark_drude_vs_standard.cpp)

# Set C++ standard
set_target_properties(benchmark_drude_vs_standard PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
)

# Enable optimization
target_compile_options(benchmark_drude_vs_standard PRIVATE -O3 -march=native)