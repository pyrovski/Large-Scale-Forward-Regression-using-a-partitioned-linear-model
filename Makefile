ifeq ($(dbg),1)
DBG=-D_DEBUG
OPT_FLAGS = -g
CUDA_FLAGS = -G
else
DBG=
OPT_FLAGS = -O3
CUDA_FLAGS = 
endif

# warning: this breaks things on other architectures...
-include CUDA_SM
CUDA_SM_VER ?= sm_13
CUDA_FLAGS+=-arch $(CUDA_SM_VER) -Xptxas -v
# -maxrregcount=16
-include CUDA_PATHS
CUDA_SDK?=/home/user/NVIDIA_GPU_Computing_SDK3.2
CUDA_TK?=/usr/local/cuda
CUDA_INC=-I$(CUDA_SDK)/C/common/inc -I$(CUDA_TK)/include

CC = $(CUDA_TK)/bin/nvcc

# The location (and name) of the BLAS/Lapack libraries
BLAS_LAPACK_LIB = -lf77blas -latlas -L/usr/lib/atlas/ -llapack
# You should not need to edit anything below this line...

# Location of GSL header files and libraries
#GSL_INCLUDE = /sw/include
GSL_LIB = -lgsl

LIBS=$(BLAS_LAPACK_LIB) $(GSL_LIB) -lcublas

# Set all include flags here for later use
#INCLUDE_FLAGS = $(GSL_INCLUDE)

SRCS = reference_glm.cu plm.cu cuda_blas.cu type.h fortran_matrix.h print_matrix.h
target = reference_glm

all: $(target)

clean:
	rm -f $(EXEC_FILES) *~ *.P *.o $(target)

$(target): $(SRCS)
	$(CC) $(DBG) $(OPT_FLAGS) $(CUDA_INC) $(CUDA_FLAGS) $< -o $@ $(LIBS)

