cmake_minimum_required(VERSION 3.12)
project(OpenMP_MPI)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(OpenMP_MPI main.cpp)

target_link_libraries(OpenMP_MPI ${MPI_C_LIBRARIES})
