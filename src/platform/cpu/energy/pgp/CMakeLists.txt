# CMakeLists.txt for PGP submodules

set(PGP_SOURCES
    PGPGrid.cpp
    PGPPrecompute.cpp
    PGPInterpolation.cpp
    PGPEvaluator.cpp
)

set(PGP_HEADERS
    PGPCore.hpp
    PGPGrid.hpp
    PGPPrecompute.hpp
    PGPInterpolation.hpp
    PGPEvaluator.hpp
)

# Add PGP sources to parent target
# These will be included in the main pygcmc_platform target
target_sources(pygcmc_platform PRIVATE ${PGP_SOURCES})

# Install headers if needed
# install(FILES ${PGP_HEADERS} DESTINATION include/platform/cpu/energy/pgp) 