#export SLEPC_DIR=/home/guido/libs/slepc-3.7.1
#export PETSC_DIR=/home/guido/libs/petsc-3.7.2
#export PETSC_ARCH=arch-linux2-c-opt

PWD:= $(shell pwd)

SRC_DIR=./src
OBJ_DIR=./obj
DEP_DIR=${PWD}/headers

CFLAGS=-g -O0
	
DEPS = ${DEP_DIR}/funct.h              \
       ${DEP_DIR}/gmsh.h               \
       ${DEP_DIR}/mesh.h               \
       ${DEP_DIR}/fem.h                \
       ${DEP_DIR}/fun.h                \
       ${DEP_DIR}/list.h               \
       ${DEP_DIR}/types.h              \
       ${DEP_DIR}/global.h             \
       ${DEP_DIR}/utils.h              

OBJ  = ${OBJ_DIR}/main.o               \
       ${OBJ_DIR}/mesh.o               \
       ${OBJ_DIR}/list.o               \
       ${OBJ_DIR}/gmsh.o               \
       ${OBJ_DIR}/fem.o                \
       ${OBJ_DIR}/fun.o                \
       ${OBJ_DIR}/ferassm.o            \
       ${OBJ_DIR}/ferinit.o            \
       ${OBJ_DIR}/ferfini.o            \
       ${OBJ_DIR}/ferboun.o            \
       ${OBJ_DIR}/ferstep.o            \
       ${OBJ_DIR}/ferelem.o            \
       ${OBJ_DIR}/fersolv.o            \
       ${OBJ_DIR}/ferpowe.o            \
       ${OBJ_DIR}/ferrods.o            \
       ${OBJ_DIR}/lst2msh.o            \
       ${OBJ_DIR}/utils.o              \
       ${OBJ_DIR}/output.o             \
       ${OBJ_DIR}/parser.o             

.PHONY: clean_
	
all: ${OBJ} 
	gcc -o fermi $^ ${SLEPC_EPS_LIB} 
	
${OBJ_DIR}/%.o: ${SRC_DIR}/%.c $(DEPS)
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I${DEP_DIR}

clean_:	    
	rm -f $(OBJ) fermi

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common



