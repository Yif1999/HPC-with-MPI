cmake_minimum_required (VERSION 2.6)
project(panel_thermal_problem)
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -g")

message("building...")

add_executable(solver
./src/2DThermal.cpp
./src/meshdef.h
./src/output_vtk.cpp
./src/solver_parallel.cpp
./src/solver_serial.cpp
)
