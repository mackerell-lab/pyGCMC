# ./CMakeLists.txt

cmake_minimum_required(VERSION 3.15)
project(pyGCMC 
    VERSION 0.1.0
    DESCRIPTION "Python-based General-purpose Classical Molecular simulation Code"
    LANGUAGES CXX
)

# Global settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type configuration
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

# Options
option(BUILD_TESTING "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(USE_SYSTEM_DEPS "Use system-installed dependencies instead of downloading them" OFF)
option(WITH_CUDA "Build with CUDA support" OFF)
option(WITH_ROCM "Build with ROCm support" OFF)
option(WITH_OPENCL "Build with OpenCL support" OFF)

# Include FetchContent
include(FetchContent)

# ---------------------------
# 1. Dependencies
# ---------------------------
# pybind11
if(USE_SYSTEM_DEPS)
    find_package(pybind11 CONFIG REQUIRED)
else()
    message(STATUS "Downloading pybind11...")
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.10.1
    )
    FetchContent_MakeAvailable(pybind11)
endif()

# Google Test (only if testing is enabled)
if(BUILD_TESTING)
    if(USE_SYSTEM_DEPS)
        find_package(GTest CONFIG REQUIRED)
    else()
        message(STATUS "Downloading Google Test...")
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.13.0
        )
        # Prevent GoogleTest from overriding our compiler/linker options
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
        set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
        
        # Ensure PIC is enabled for GTest
        set_target_properties(gtest gtest_main PROPERTIES 
            POSITION_INDEPENDENT_CODE ON
        )
    endif()
endif()

# CUDA配置
if(WITH_CUDA)
    enable_language(CUDA)
    find_package(CUDA REQUIRED)
    if(CUDA_FOUND)
        message(STATUS "CUDA found: ${CUDA_VERSION}")
        add_definitions(-DWITH_CUDA)
    else()
        message(FATAL_ERROR "CUDA requested but not found.")
    endif()
endif()

# ROCm配置
if(WITH_ROCM)
    find_package(hip REQUIRED)
    if(hip_FOUND)
        message(STATUS "ROCm found: ${hip_VERSION}")
        add_definitions(-DWITH_ROCM)
    else()
        message(FATAL_ERROR "ROCm requested but not found.")
    endif()
endif()

# OpenCL配置
if(WITH_OPENCL)
    find_package(OpenCL REQUIRED)
    if(OpenCL_FOUND)
        message(STATUS "OpenCL found: ${OpenCL_VERSION}")
        add_definitions(-DWITH_OPENCL)
    else()
        message(FATAL_ERROR "OpenCL requested but not found.")
    endif()
endif()

# ---------------------------
# 2. Core Library
# ---------------------------
add_subdirectory(modules/core)

# ---------------------------
# 3. GCMC Module
# ---------------------------
# add_subdirectory(modules/gcmc)

# ---------------------------
# 4. Optimization Module
# ---------------------------
# add_subdirectory(modules/optimization)

# ---------------------------
# 5. Bindings
# ---------------------------
add_subdirectory(modules/bindings)

# ---------------------------
# 6. Tests
# ---------------------------
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---------------------------
# 7. Examples
# ---------------------------
if(BUILD_EXAMPLES)
    add_subdirectory(examples/basic)
    # add_subdirectory(examples/advanced)
endif()

# ---------------------------
# 8. Installation
# ---------------------------
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/pyGCMC)

install(TARGETS pygcmc_core
    EXPORT pyGCMC-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT pyGCMC-targets
    FILE pyGCMCTargets.cmake
    NAMESPACE pyGCMC::
    DESTINATION ${INSTALL_CONFIGDIR}
)

# Create version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/pyGCMCConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install version file
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/pyGCMCConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
