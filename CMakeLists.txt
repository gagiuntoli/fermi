cmake_minimum_required(VERSION 3.0)

set(CMAKE_C_COMPILER "/usr/bin/mpicc")
project(Fermi VERSION 0.1 LANGUAGES C)

if(EXISTS "$ENV{PETSC_DIR}/include/" AND EXISTS "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include/")
	set(PETSC_INC "$ENV{PETSC_DIR}/include" "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include")
else()
	message(SEND_ERROR "PETSc includes not found")
endif()

if(EXISTS "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/libpetsc.so")
	set(PETSC_LIB "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/libpetsc.so")
else()
	message(SEND_ERROR "PETSc library not found")
endif()

find_package(MPI REQUIRED)

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR} ${PETSC_INC})

# wildcard all the sources in src
file(GLOB SOURCES src/*.c)

link_libraries("-lm" ${PETSC_LIB})
#target_link_libraries(${SLEPC_LIB})

# Executable
add_executable(fermi ${SOURCES})
