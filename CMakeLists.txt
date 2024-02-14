cmake_minimum_required(VERSION 3.20)

project(Fermi VERSION 0.1 LANGUAGES C)

find_package(MPI REQUIRED)

if(NOT EXISTS "$ENV{PETSC_DIR}/include/" OR NOT EXISTS "$ENV{PETSC_DIR}/lib/")
 	message(SEND_ERROR "PETSc not found")
endif()
message("MPI " ${MPI_INCLUDE_PATH})

include_directories(include ${MPI_INCLUDE_PATH} $ENV{PETSC_DIR}/include)

link_libraries("-lm" $ENV{PETSC_DIR}/lib/libpetsc_real.so ${MPI_C_LIBRARIES})

set(SOURCES
	src/gmsh.c 
	src/list.c 
	src/solve.c 
	src/utils.c 
	src/rods.c 
	src/boundary.c 
	src/finish.c 
	src/assembly.c 
	src/elemental.c 
	src/fem.c 
	src/finish.c 
	src/init.c 
	src/fun.c 
	src/lst2msh.c 
	src/mesh.c 
	src/output.c 
	src/parser.c 
	src/power.c 
)

add_library(fermi_lib SHARED ${SOURCES})

add_executable(fermi_sim src/main.c)

target_link_libraries(fermi_sim fermi_lib)

# Tests
add_executable(test_parser src/test_parser.c)
target_link_libraries(test_parser fermi_lib)
add_test(NAME test_1 COMMAND test_parser)

enable_testing()
