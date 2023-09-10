#################################################################
#  CMake configure file for Monte Carlo eXtreme (mmc)
#  Qianqian Fang <q.fang at neu.edu>
#  2018/08/26
#################################################################
# Note: when multiple gcc versions are installed, cmake may find
# the highest version, but nvcc may throw an error saying it does
# not support this gcc version. To solve this, one should run
#      CC=gcc-X CXX=g++-X cmake ..
# where "X" is the version (such as 9) that both installed on
# your system and is supported by nvcc
#################################################################

cmake_minimum_required(VERSION 3.3)

project(mmc)

if(${CMAKE_VERSION} VERSION_GREATER "3.11.20")
    cmake_policy(SET CMP0074 NEW)  # to make -DOPENCL_ROOT= work with newer CMake versions as well
endif()

find_package(OpenMP REQUIRED)
find_package(OpenCL REQUIRED)

add_subdirectory(zmat)

option(BUILD_MEX "Build mex" ON)
option(BUILD_CUDA "Build cuda" OFF)

if(BUILD_PYTHON)
    add_subdirectory(pybind11)
    find_package (Python3 COMPONENTS Interpreter Development)
    include_directories(${PYTHON_INCLUDE_DIRS})
endif()

if(BUILD_MEX)
    find_package(Matlab)
endif()

if(BUILD_CUDA)
  find_package(CUDA QUIET REQUIRED)
endif()

# C Options
set(CMAKE_CXX_FLAGS "-Wall -g -DMCX_EMBED_CL -fno-strict-aliasing -m64 -DMMC_USE_SSE -DHAVE_SSE2 -msse -msse2 -msse3 -mssse3 -msse4.1 -O3 -DUSE_OS_TIMER -DUSE_OPENCL -DMMC_XORSHIFT -fPIC ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(RESOURCE_COMPILER xxd)
set(RESOURCE_INPUT "mmc_core.cl")
set(RESOURCE_OUTPUT "${CMAKE_SOURCE_DIR}/mmc_core.clh")

add_custom_target(
    clheader
    COMMAND bash -c "cd ${CMAKE_SOURCE_DIR}; ${RESOURCE_COMPILER} -i ${RESOURCE_INPUT} | sed 's/\\([0-9a-f]\\)$/\\0, 0x00/' > ${RESOURCE_OUTPUT}"
    comment "Compiling ${RESOURCE_INPUT} to ${RESOURCE_OUTPUT}"
    VERBATIM
)

if(BUILD_CUDA)
  # NVC# NVCC Options
  set(
      CUDA_NVCC_FLAGS
      ${CUDA_NVCC_FLAGS};
      -g -lineinfo -DUSE_ATOMIC -DMCX_SAVE_DETECTORS -DMCX_DO_REFLECTION -DUSE_DMMC -DUSE_BLBADOUEL -Xcompiler -fPIC -use_fast_math -Xptxas -O3,-v -Xcompiler -Wall -Xcompiler "${OpenMP_CXX_FLAGS}"
      )
endif()

# Add include directories
include_directories(cjson ubj zmat zmat/easylzma)

# Add link directories
link_directories(zmat)

# Create mmc library
add_library(mmc STATIC
    mmc_utils.c
    mmc_utils.h
    mmc_mesh.c
    mmc_mesh.h
    mmc_host.c
    mmc_host.h
    mmc_raytrace.c
    mmc_raytrace.h
    mmc_rand_xorshift128p.c
    mmc_rand_xorshift128p.h
    mmc_bench.h
    mmc_tictoc.c
    mmc_tictoc.h
    mmc_cl_utils.c
    mmc_cl_utils.h
    mmc_cl_host.c
    mmc_cl_host.h
    mmc_highorder.cpp
    cjson/cJSON.c
    cjson/cJSON.h
    ubj/ubj.h
    ubj/ubjw.c
    )

# Add all project units
add_executable(
    mmc-exe
    mmc.c
    )

set_target_properties(mmc-exe
        PROPERTIES OUTPUT_NAME mmc)

# Link options
target_link_libraries(
    mmc-exe
    mmc OpenMP::OpenMP_CXX
    zmat
    m
    OpenCL::OpenCL
    )

add_dependencies(mmc clheader)
add_dependencies(mmc-exe clheader)

#static link libraries
set_target_properties(mmc-exe PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(mmc-exe PROPERTIES LINK_SEARCH_END_STATIC 1)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
endif()

if(BUILD_CUDA)
    # Create mmc library
    cuda_add_library(mmc-cuda STATIC
        mmc_cu_host.cu
        mmc_cu_host.h
        mmc_utils.c
        mmc_utils.h
        )

  # Create mmc-trinity binary
  cuda_add_executable(
      mmc-trinity
      mmc.c
      )
  # Link options
  target_link_libraries(
      mmc-trinity
      mmc mmc-cuda OpenMP::OpenMP_CXX
      zmat
      m
      OpenCL::OpenCL
      )

  target_compile_definitions(mmc-trinity PUBLIC USE_CUDA)
  add_dependencies(mmc-trinity mmc mmc-cuda)

  set_target_properties(mmc-trinity
        PROPERTIES OUTPUT_NAME mmciii)
  #static link libraries
  set_target_properties(mmc-trinity PROPERTIES LINK_SEARCH_START_STATIC 1)
  set_target_properties(mmc-trinity PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

# Build mex file
if(BUILD_MEX AND Matlab_FOUND)
    # Create mcx-matlab library
    add_library(mmc-matlab STATIC
            mmc_utils.c
            mmc_utils.h
            mmc_mesh.c
            mmc_mesh.h
            mmc_host.c
            mmc_host.h
            mmc_raytrace.c
            mmc_raytrace.h
            mmc_rand_xorshift128p.c
            mmc_rand_xorshift128p.h
            mmc_bench.h
            mmc_tictoc.c
            mmc_tictoc.h
            mmc_cl_utils.c
            mmc_cl_utils.h
            mmc_cl_host.c
            mmc_cl_host.h
            mmc_highorder.cpp
            cjson/cJSON.c
            cjson/cJSON.h
            ubj/ubj.h
            ubj/ubjw.c
            )

    target_compile_definitions(mmc-matlab PUBLIC MCX_CONTAINER MATLAB_MEX_FILE)

    if(${CMAKE_VERSION} VERSION_LESS "3.24.0")
            matlab_add_mex(
              NAME mmclab
              SRC mmclab.cpp
              LINK_TO mmc-matlab OpenMP::OpenMP_CXX OpenCL::OpenCL zmat
            )
    else()
            matlab_add_mex(
              NAME mmclab
              SRC mmclab.cpp
              NO_IMPLICIT_LINK_TO_MATLAB_LIBRARIES
              LINK_TO ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} mmc-matlab OpenMP::OpenMP_CXX OpenCL::OpenCL zmat
            )
    endif()

    target_compile_definitions(mmclab PUBLIC MCX_CONTAINER MATLAB_MEX_FILE)

    set_target_properties(mmclab
            PROPERTIES OUTPUT_NAME ${CMAKE_SOURCE_DIR}/../mmclab/mmc)
    add_dependencies(mmclab clheader)

    if(BUILD_CUDA)
        if(${CMAKE_VERSION} VERSION_LESS "3.24.0")
                matlab_add_mex(
                  NAME mmclab-trinity
                  SRC mmclab.cpp
                  LINK_TO mmc-cuda mmc-matlab OpenMP::OpenMP_CXX OpenCL::OpenCL zmat
                )
        else()
                matlab_add_mex(
                  NAME mmclab-trinity
                  SRC mmclab.cpp
                  NO_IMPLICIT_LINK_TO_MATLAB_LIBRARIES
                  LINK_TO ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} mmc-cuda mmc-matlab OpenMP::OpenMP_CXX OpenCL::OpenCL zmat
                )
        endif()

        target_compile_definitions(mmclab-trinity PUBLIC MCX_CONTAINER MATLAB_MEX_FILE USE_CUDA)

        set_target_properties(mmclab-trinity
                PROPERTIES OUTPUT_NAME ${CMAKE_SOURCE_DIR}/../mmclab/mmciii)
        add_dependencies(mmclab-trinity clheader)
    endif()


endif()
