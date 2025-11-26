cmake_minimum_required(VERSION 3.10)
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER g++)
endif()

execute_process(
  COMMAND git log --pretty=format:%h -n 1
  OUTPUT_VARIABLE GIT_COMMIT
)
if ("${GIT_COMMIT}" STREQUAL "")
  set(GIT_COMMIT "NA")
endif()

file(STRINGS "VERSION" RM_VERSION)
project(remeta VERSION ${RM_VERSION})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -Wall -ffast-math -fvisibility=hidden -fopenmp")
add_definitions(-DVERSION_NUMBER="v${RM_VERSION}-${GIT_COMMIT}")

add_executable(remeta
  ${CMAKE_SOURCE_DIR}/src/io/allele_freq_writer.cpp
  ${CMAKE_SOURCE_DIR}/src/io/anno_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/bgz_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/bgz_writer.cpp
  ${CMAKE_SOURCE_DIR}/src/io/block_pgen_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/gene_set_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/htpv4_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/htpv4_writer.cpp
  ${CMAKE_SOURCE_DIR}/src/io/remeta_matrix_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/regenie_ld_matrix_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/regenie_ld_matrix_writer.cpp
  ${CMAKE_SOURCE_DIR}/src/io/ref_ld_block_mapper.cpp
  ${CMAKE_SOURCE_DIR}/src/io/ref_ld_genetic_map.cpp
  ${CMAKE_SOURCE_DIR}/src/io/ref_ld_matrix_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/ref_ld_matrix_writer.cpp
  ${CMAKE_SOURCE_DIR}/src/io/regenie_anno_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/io/tabixed_anno_reader.cpp
  ${CMAKE_SOURCE_DIR}/src/meta/es_meta_analyzer.cpp
  ${CMAKE_SOURCE_DIR}/src/meta/genep_meta_analyzer.cpp
  ${CMAKE_SOURCE_DIR}/src/meta/htp_meta_analyzer.cpp
  ${CMAKE_SOURCE_DIR}/src/meta/pv_meta_analyzer.cpp
  ${CMAKE_SOURCE_DIR}/src/meta/set_meta_analyzer.cpp
  ${CMAKE_SOURCE_DIR}/src/meta/variant_meta_analyzer.cpp
  ${CMAKE_SOURCE_DIR}/src/stat/misc.cpp
  ${CMAKE_SOURCE_DIR}/src/stat/spa.cpp
  ${CMAKE_SOURCE_DIR}/src/stat/tests.cpp
  ${CMAKE_SOURCE_DIR}/src/allele_frequency_map.cpp
  ${CMAKE_SOURCE_DIR}/src/annotation_map.cpp
  ${CMAKE_SOURCE_DIR}/src/htpv4_pos.cpp
  ${CMAKE_SOURCE_DIR}/src/logging.cpp
  ${CMAKE_SOURCE_DIR}/src/parameter_checks.cpp
  ${CMAKE_SOURCE_DIR}/src/remeta.cpp
  ${CMAKE_SOURCE_DIR}/src/run_compute_ref_ld.cpp
  ${CMAKE_SOURCE_DIR}/src/run_esma.cpp
  ${CMAKE_SOURCE_DIR}/src/run_genep.cpp
  ${CMAKE_SOURCE_DIR}/src/run_htp.cpp
  ${CMAKE_SOURCE_DIR}/src/run_index_anno.cpp
  ${CMAKE_SOURCE_DIR}/src/run_ld_deflate.cpp
  ${CMAKE_SOURCE_DIR}/src/run_ld_inflate.cpp
  ${CMAKE_SOURCE_DIR}/src/run_pvma.cpp
  ${CMAKE_SOURCE_DIR}/src/source_map.cpp
  ${CMAKE_SOURCE_DIR}/src/util.cpp
  ${CMAKE_SOURCE_DIR}/src/variant_filter.cpp
)

find_path(
  EIGEN
  Eigen
  HINTS ${EIGEN_PATH}
        /usr/local/include/eigen3
        /usr/local/lib/eigen-3.4.0
)
if (${EIGEN} STREQUAL EIGEN-NOTFOUND)
  message(SEND_ERROR "could not find Eigen folder at EIGEN_PATH")
else()
  message("EIGEN: ${EIGEN}")
endif()

target_include_directories(remeta PRIVATE ${EIGEN})

find_library(
  HTSLIB
  libhts.a
  HINTS ${HTSLIB_PATH}
)
find_library(
  DEFLATE
  libdeflate.a
)
if (${HTSLIB} STREQUAL HTSLIB-NOTFOUND)
  message(SEND_ERROR "Could not find libhts.a at ${HTSLIB_PATH}")
else()
  message("HTSlib: ${HTSLIB}")
endif()
if (${DEFLATE} STREQUAL DEFLATE-NOTFOUND)
  message(SEND_ERROR "Could not find libdeflate.a")
endif()

find_library(
  BOOST_LOG_SETUP
  libboost_log_setup.a
  HINTS ${BOOST_PATH}
)
find_library(
  BOOST_LOG
  libboost_log.a
  HINTS ${BOOST_PATH}
)
find_library(
  BOOST_PROGRAM_OPTIONS
  libboost_program_options.a
  HINTS ${BOOST_PATH}
)
find_library(
  BOOST_FILESYSTEM
  libboost_filesystem.a
  HINTS ${BOOST_PATH}
)
find_library(
  BOOST_SYSTEM
  libboost_system.a
  HINTS ${BOOST_PATH}
)
find_library(
  BOOST_THREAD
  libboost_thread.a
  HINTS ${BOOST_PATH}
)
find_library(Z libz.so)
find_library(SSL libssl.so)
find_library(CRYPTO libcrypto.so)
find_library(BZ2 libbz2.so)
find_library(LZMA liblzma.so)
find_library(CURL libcurl.so)

set(EXTERNAL_LIBS
  ${HTSLIB}
  ${DEFLATE}
  ${BOOST_LOG_SETUP}
  ${BOOST_LOG}
  ${BOOST_PROGRAM_OPTIONS}
  ${BOOST_FILESYSTEM}
  ${BOOST_SYSTEM}
  ${BOOST_THREAD}
  ${PTHREAD}
  ${Z}
  ${SSL}
  ${CRYPTO}
  ${DL}
  ${BZ2}
  ${LZMA}
  ${CURL}
)

target_link_libraries(remeta PUBLIC ${EXTERNAL_LIBS})

add_custom_target(
  pgenlib.a
  COMMAND make
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/pgenlib
)
target_link_libraries(
  remeta PUBLIC
  ${CMAKE_SOURCE_DIR}/lib/pgenlib/pgenlib.a
)
target_include_directories(
  remeta PUBLIC
  ${CMAKE_SOURCE_DIR}/lib/pgenlib
)
add_dependencies(remeta pgenlib.a)

add_custom_target(
  faddeeva.a
  COMMAND make
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/faddeeva
)
target_link_libraries(
  remeta PUBLIC
  ${CMAKE_SOURCE_DIR}/lib/faddeeva/faddeeva.a
)
target_include_directories(
  remeta PUBLIC
  ${CMAKE_SOURCE_DIR}/lib/faddeeva
)
add_dependencies(remeta faddeeva.a)

add_custom_target(
  qfc.a
  COMMAND make
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/qfc
)
target_link_libraries(
  remeta PUBLIC
  ${CMAKE_SOURCE_DIR}/lib/qfc/qfc.a
)
target_include_directories(
  remeta PUBLIC
  ${CMAKE_SOURCE_DIR}/lib/qfc
)
add_dependencies(remeta qfc.a)

#####################################
### Link Linear Algebra Libraries ###
#####################################
# Intel MKL
find_library(MKL libmkl_intel_lp64.a
  HINTS ${MKLROOT} $ENV{MKLROOT}/lib
)
find_library(MKL_GNU_THREAD libmkl_gnu_thread.a
  HINTS ${MKLROOT} $ENV{MKLROOT}/lib
)
find_library(MKL_CORE libmkl_core.a
  HINTS ${MKLROOT} $ENV{MKLROOT}/lib
)
find_library(GOMP libgomp.a HINTS /usr/lib/gcc/x86_64-linux-gnu/8/)
find_library(M NAMES m)

# Openblas
find_library(LAPACK liblapack.a)
find_library(LAPACKE liblapacke.a)
find_library(OPENBLAS libopenblas.a)

if (MKL AND MKL_GNU_THREAD AND MKL_CORE)
  message("Compiling with Intel MKL")
  target_include_directories(remeta PRIVATE ${MKL_INCLUDE})
  message("\t${MKL}")
  message("\t${MKL_GNU_THREAD}")
  message("\t${MKL_CORE}")
  target_link_libraries(remeta PRIVATE
    "-Wl,--start-group"
    ${MKL} ${MKL_GNU_THREAD} ${MKL_CORE}
    "-Wl,--end-group"
    -lgomp
    ${M}
    -ldl)
  #set(MKL_LINK static)
  find_package(MKL CONFIG PATHS ${MKLROOT} $ENV{MKLROOT}} QUIET)
  target_include_directories(remeta PRIVATE ${MKL_INCLUDE})
  target_compile_definitions(remeta PUBLIC "-DMKL_ILP64 -DWITH_MKL -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE")
elseif (LAPACK AND LAPACKE AND OPENBLAS)
  message("Compiling with OpenBLAS")
  message("\t${LAPACK}")
  message("\t${LAPACKE}")
  message("\t${OPENBLAS}")
  set(BLAS_LIBS ${LAPACK} ${LAPACKE} ${OPENBLAS})

  target_link_libraries(remeta PUBLIC ${BLAS_LIBS})
  target_compile_definitions(remeta PUBLIC "-DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE")
else()
  message("Missing libraries for Intel MKL and OpenBLAS")
  message("\tLAPACK ${LAPACK}")
  message("\tLAPACKE ${LAPACKE}")
  message("\tOPENBLAS ${OPENBLAS}")
  message("\tMKL ${MKL}")
  message("\tMKL_GNU_THREAD ${MKL_GNU_THREAD}")
  message("\tMKL_CORE ${MKL_CORE}")
endif()