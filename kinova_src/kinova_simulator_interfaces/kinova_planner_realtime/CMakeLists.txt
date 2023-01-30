project(PZsparse-pybind)

cmake_minimum_required(VERSION 3.4...3.18)
find_package(CUDA REQUIRED)
find_package(PythonInterp 3.6 REQUIRED)
find_package(PythonLibs 3.6 REQUIRED)
find_package(OpenMP REQUIRED)

set(IPOPT_LIBRARY_DIRS /usr/local/lib)
set(IPOPT_INCLUDE_DIRS /usr/local/include/coin-or)

include_directories(
    ${IPOPT_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    )

link_directories(
    ${IPOPT_LIBRARY_DIRS}
    /usr/local/cuda/lib64
    )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++14")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_subdirectory(pybind11)

cuda_add_library(rtd_force SHARED
  rtd_force_pybind.cu
  NLPclass.cu
  # NLPclass_freeinput.cu
  CollisionChecking.cu
  Dynamics.cu
  Trajectory.cu
  PZsparse.cu
  )

target_link_libraries(rtd_force # PUBLIC pybind11::module
  ${PYTHON_LIBRARIES}
  ipopt
  cuda
  cudart
  OpenMP::OpenMP_CXX
  )

# pybind11_add_module(armtd_main_pybind SHARED
#   armtd_main_pybind.cpp
#   NLPclass.cu
#   CollisionChecking.cu
#   Dynamics.cpp
#   Trajectory.cpp
#   PZsparse.cpp
#   )

pybind11_extension(rtd_force)
# set_target_properties(armtd_main_pybind PROPERTIES PREFIX "")
# set_target_properties(test_pybind PROPERTIES PREFIX "")

# if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
#     pybind11_strip(example)
# endif()

# set_target_properties(armtd_main_pybind PROPERTIES CXX_VISIBILITY_PRESET "hidden"
#                                          CUDA_VISIBILITY_PRESET "hidden")
