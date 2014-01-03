# cmake version
cmake_minimum_required (VERSION 2.8)

# project name
project (LinBolt)

# path to FindGSL.cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

# require GSL
find_package(GSL REQUIRED)

# build
add_executable(LinBolt arsenal.cpp medium.cpp particle.cpp system.cpp main.cpp ParameterReader.cpp particle.cpp scattering.cpp routines.cpp checks.cpp)

target_link_libraries(LinBolt ${GSL_LIBRARIES}  ${GSLCBLAS_LIBRARIES})
