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
CUDA_LIBS=-L/usr/local/cuda/lib64 -lcudart

GPUCC = $(CUDA_TK)/bin/nvcc
CC=mpicxx

# The location (and name) of the BLAS/Lapack libraries
BLAS_LAPACK_LIB = -llapack
# You should not need to edit anything below this line...

# Location of GSL header files and libraries
#GSL_INCLUDE = /sw/include
GSL_LIB = -lgsl

LIBS=$(BLAS_LAPACK_LIB) $(GSL_LIB) -lcublas $(CUDA_LIBS)

# Set all include flags here for later use
#INCLUDE_FLAGS = $(GSL_INCLUDE)

CPU_SRC = reference_glm.cpp tvUtil.cpp
GPU_SRC = plm.cu cuda_blas.cu
HEADERS= type.h fortran_matrix.h print_matrix.h glm.h plm.h
SRC=$(CPU_SRC) $(GPU_SRC)
target = reference_glm
objects = $(patsubst %.cpp,%.o,$(CPU_SRC)) $(patsubst %.cu,%.o,$(GPU_SRC))
all: $(target)

clean:
	rm -f $(EXEC_FILES) *~ *.P *.o $(target)

$(target): $(objects)
#	$(CC) $(DBG) $(OPT_FLAGS) $(CUDA_INC) $(CPU_SRC) -c 
#	$(GPUCC) $(DBG) $(OPT_FLAGS) $(CUDA_INC) $(GPU_SRC) $(CUDA_FLAGS) $(LIBS)
	$(CC) -o $@ $^ $(LIBS)

%.o:%.cpp
	$(CC) $(DBG) $(OPT_FLAGS) $(CUDA_INC) $^ -c
%.o:%.cu
	$(GPUCC) $(DBG) $(OPT_FLAGS) $(CUDA_INC) $(CUDA_FLAGS) $^ -c