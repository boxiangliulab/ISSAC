cmake_minimum_required(VERSION 3.10)

project(MyProject VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(gtf/
                    utils/
                    utils/htslib/
                    utils/htslib/htslib/
                    utils/bedtools/bedFile/
                    utils/bedtools/lineFileUtilities/
                    utils/bedtools/gzstream/
                    utils/bedtools/fileType/
                    utils/bedtools/stringUtilities/
                    utils/samtools/)

# Find the HTSlib package installed by Conda
find_path(HTSLIB_INCLUDE_DIR htslib/hts.h HINTS $ENV{CONDA_PREFIX}/include)
find_library(HTSLIB_LIBRARY NAMES hts HINTS $ENV{CONDA_PREFIX}/lib)
find_library(CRYPTO_LIB NAMES crypto HINTS $ENV{CONDA_PREFIX}/lib)

find_path(GSL_INCLUDE_DIR gsl HINTS $ENV{CONDA_PREFIX}/include)
find_library(GSL_LIBRARY NAMES gsl HINTS $ENV{CONDA_PREFIX}/lib)

set(NLOPT_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")
set(NLOPT_LIBRARY_DIR "$ENV{CONDA_PREFIX}/lib/libnlopt.so.0.12.0")
# Handle the case where HTSlib was not found
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
    message(FATAL_ERROR "HTSlib not found in your Conda environment.")
endif()

find_path(EIGEN_INCLUDE_DIR eigen3 HINTS $ENV{CONDA_PREFIX}/include)
find_library(EIGEN_LIBRARY NAMES eigen HINTS $ENV{CONDA_PREFIX}/lib)

find_package(Eigen3 REQUIRED)
# Include HTSlib headers and link the library
include_directories(${HTSLIB_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${NLOPT_INCLUDE_DIR})
#The main executable
add_executable (ISSAC eigen_decom.cpp eigen_decom.h IR_denominator_quant.cpp pheno.h phenotype_group.cpp phenotype_output.cpp 
  IR_denominator_quant.h junc_stat.cpp junctions_extractor.cc junctions_extractor.h main.cpp)

target_link_libraries(ISSAC ${HTSLIB_LIBRARY})

target_link_libraries(ISSAC ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARIES})

target_link_libraries(ISSAC Eigen3::Eigen)

target_link_libraries(ISSAC ${NLOPT_LIBRARY_DIR})

target_link_libraries(ISSAC ${CRYPTO_LIB})
