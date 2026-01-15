##############################################################################
#
#  unix-mpicc-default.mk
#
#  Compile for parallel execution assuming an mpi compiler
#  wrapper "mpicc" is available.
#
##############################################################################

BUILD   = parallel
MODEL   = -D_D3Q19_
TARGET  =

CC      = mpicc -fopenmp
CFLAGS  = -O2 -g -Wall 

AR      = ar                      # standard ar command
ARFLAGS = -cru                    # flags for ar
LDFLAGS =                         # additional link time flags

MPI_INC_PATH      =  -I$(PETSC_DIR)/$(PETSC_ARCH)/include/             # path to mpi.h (if required)
MPI_LIB_PATH      =  -L$(PETSC_DIR)/$(PETSC_ARCH)/lib/           # path to libmpi.a (if required)
MPI_LIB           =  -lmpi            # -lmpi (if required)

LAUNCH_MPIRUN_CMD = mpirun -np 1

HAVE_PETSC	=	true
PETSC_INC = -I$(PETSC_DIR)/include/ -I$(PETSC_DIR)/$(PETSC_ARCH)/include/   
PETSC_LIB =   -lpetsc #