cmake_minimum_required(VERSION 3.20)

project(Fermi VERSION 0.1 LANGUAGES CXX)

find_package(MPI REQUIRED)

if(NOT EXISTS "$ENV{PETSC_DIR}/include/" OR NOT EXISTS "$ENV{PETSC_DIR}/lib/")
 	message(SEND_ERROR "PETSc not found")
endif()
message("MPI " ${MPI_INCLUDE_PATH})

include_directories(include ${MPI_INCLUDE_PATH} $ENV{PETSC_DIR}/include submodules/tomlplusplus)

link_libraries("-lm" $ENV{PETSC_DIR}/lib/libpetsc_real.so ${MPI_CXX_LIBRARIES})

set(SOURCES
	src/gmsh.cpp 
	src/list.cpp 
	src/solve.cpp 
	src/utils.cpp 
	src/rods.cpp 
	src/boundary.cpp 
	src/finish.cpp 
	src/assembly.cpp 
	src/elemental.cpp 
	src/fem.cpp 
	src/finish.cpp 
	src/init.cpp 
	src/fun.cpp 
	src/lst2msh.cpp 
	src/mesh.cpp 
	src/output.cpp 
	src/parser.cpp 
	src/power.cpp 
)


add_library(fermi_lib SHARED ${SOURCES})

add_executable(fermi_sim src/main.cpp)

target_link_libraries(fermi_sim PUBLIC fermi_lib)

# Tests
add_executable(test_parser src/test_parser.cpp)
target_link_libraries(test_parser fermi_lib)
add_test(NAME test_1 COMMAND test_parser)

enable_testing()
