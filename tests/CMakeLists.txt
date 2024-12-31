# tests/CMakeLists.txt

# Enable testing
enable_testing()

# Set test data directory
set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

# Define paths to parameter files
set(WATER_PRM_FILE "${TEST_DATA_DIR}/toppar_water_ions.str")
set(PROT_PRM_FILE "${TEST_DATA_DIR}/par_all36m_prot.prm")

# Add compile definitions
add_compile_definitions(TEST_DATA_DIR="${TEST_DATA_DIR}")

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/modules/core/include
    ${CMAKE_BINARY_DIR}/tests
)

# Add Test Executables
add_executable(test_pdb_parser core/test_pdb_parser.cpp)
target_link_libraries(test_pdb_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_top_parser core/test_top_parser.cpp)
target_link_libraries(test_top_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_psf_parser core/test_psf_parser.cpp)
target_link_libraries(test_psf_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_itp_parser core/test_itp_parser.cpp)
target_link_libraries(test_itp_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_top_and_itp core/test_top_and_itp_parser.cpp)
target_link_libraries(test_top_and_itp PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_psf_and_itp core/test_psf_and_itp_parser.cpp)
target_link_libraries(test_psf_and_itp PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_ff_parser core/test_ff_parser.cpp)
target_link_libraries(test_ff_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_pdb_top_ff_parser core/test_pdb_top_ff_parser.cpp)
target_link_libraries(test_pdb_top_ff_parser PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

add_executable(test_system core/test_system.cpp)
target_link_libraries(test_system PRIVATE 
    pygcmc_core 
    GTest::gtest_main
)

# Register Tests with CTest
add_test(NAME test_pdb_parser COMMAND test_pdb_parser)
add_test(NAME test_top_parser COMMAND test_top_parser)
add_test(NAME test_psf_parser COMMAND test_psf_parser)
add_test(NAME test_itp_parser COMMAND test_itp_parser)
add_test(NAME test_top_and_itp COMMAND test_top_and_itp)
add_test(NAME test_psf_and_itp COMMAND test_psf_and_itp)

# Register test_ff_parser with required arguments
add_test(
    NAME test_ff_parser
    COMMAND test_ff_parser
        --prm_file=${WATER_PRM_FILE}
        --prot_prm_file=${PROT_PRM_FILE}
)

add_test(NAME test_pdb_top_ff_parser COMMAND test_pdb_top_ff_parser)
add_test(NAME test_system COMMAND test_system)

# Configure config.hpp
configure_file(
    "${CMAKE_SOURCE_DIR}/tests/config.hpp.in"
    "${CMAKE_BINARY_DIR}/tests/config.hpp"
    @ONLY
)

# Add include directories for all test targets
target_include_directories(test_pdb_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_top_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_psf_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_itp_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_top_and_itp PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_psf_and_itp PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_ff_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_pdb_top_ff_parser PRIVATE ${CMAKE_BINARY_DIR}/tests)
target_include_directories(test_system PRIVATE ${CMAKE_BINARY_DIR}/tests)

# Copy parameter files to build directory
if(EXISTS "${WATER_PRM_FILE}")
    file(COPY "${WATER_PRM_FILE}" DESTINATION "${CMAKE_BINARY_DIR}/tests/data")
else()
    message(FATAL_ERROR "Parameter file toppar_water_ions.str not found in tests/data")
endif()

if(EXISTS "${PROT_PRM_FILE}")
    file(COPY "${PROT_PRM_FILE}" DESTINATION "${CMAKE_BINARY_DIR}/tests/data")
else()
    message(FATAL_ERROR "Parameter file par_all36m_prot.prm not found in tests/data")
endif()
