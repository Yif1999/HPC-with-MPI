cmake_minimum_required (VERSION 2.6)
project(panel_thermal_problem)
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")

message("building...")

add_executable(solver
./src/miniCFD.cpp
./src/util.cpp
./src/util.h
)
