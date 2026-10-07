cmake_minimum_required(VERSION 3.10)

project(ISSAC VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Project-bundled headers (the bundled htslib headers are intentionally NOT included,
# so that the htslib headers always match the linked conda library)
include_directories(gtf/
                    utils/
                    utils/bedtools/bedFile/
                    utils/bedtools/lineFileUtilities/
                    utils/bedtools/gzstream/
                    utils/bedtools/fileType/
                    utils/bedtools/stringUtilities/
                    utils/samtools/)

# ---------------- Dependencies (from the active conda environment) ----------------
set(DEP_HINT $ENV{CONDA_PREFIX})

# HTSlib
find_path(HTSLIB_INCLUDE_DIR htslib/hts.h HINTS ${DEP_HINT}/include)
find_library(HTSLIB_LIBRARY NAMES hts HINTS ${DEP_HINT}/lib)

# GSL (gsl + gslcblas)
find_path(GSL_INCLUDE_DIR gsl/gsl_cdf.h HINTS ${DEP_HINT}/include)
find_library(GSL_LIBRARY NAMES gsl HINTS ${DEP_HINT}/lib)
find_library(GSL_CBLAS_LIBRARY NAMES gslcblas HINTS ${DEP_HINT}/lib)

# NLopt (found dynamically; no hard-coded version so it works on amd64 and arm64)
find_path(NLOPT_INCLUDE_DIR nlopt.h HINTS ${DEP_HINT}/include)
find_library(NLOPT_LIBRARY NAMES nlopt HINTS ${DEP_HINT}/lib)

# OpenSSL crypto
find_library(CRYPTO_LIB NAMES crypto HINTS ${DEP_HINT}/lib)

# Eigen (header-only)
find_package(Eigen3 REQUIRED)

# Fail early with a clear message if anything is missing
foreach(dep HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY GSL_INCLUDE_DIR GSL_LIBRARY GSL_CBLAS_LIBRARY
            NLOPT_INCLUDE_DIR NLOPT_LIBRARY CRYPTO_LIB)
    if(NOT ${dep})
        message(FATAL_ERROR "${dep} not found. Please activate the conda environment containing the dependencies.")
    endif()
endforeach()

include_directories(${HTSLIB_INCLUDE_DIR} ${GSL_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR})

# ---------------- Main executable ----------------
add_executable(ISSAC
    eigen_decom.cpp eigen_decom.h
    IR_denominator_quant.cpp IR_denominator_quant.h
    pheno.h phenotype_group.cpp phenotype_output.cpp
    IR_combine.cpp junc_stat.cpp
    junctions_extractor.cc junctions_extractor.h
    main.cpp)

target_link_libraries(ISSAC
    ${HTSLIB_LIBRARY}
    ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}
    Eigen3::Eigen
    ${NLOPT_LIBRARY}
    ${CRYPTO_LIB})
