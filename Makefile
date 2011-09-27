ifeq ($(dbg),1)
DBG=-D_DEBUG
OPT_FLAGS = -g
CUDA_FLAGS = -G
else
ifeq ($(datafiles),1)
DBG=-D_DEBUG
else
DBG=
endif
OPT_FLAGS = -O3
CUDA_FLAGS = 
endif

-include CUDA_SM
-include LIBRARY_PATHS

# warning: this breaks things on other architectures...
CUDA_SM_VER ?= sm_13

CUDA_FLAGS+=-arch $(CUDA_SM_VER) -Xptxas -v
# -maxrregcount=16
-include CUDA_PATHS
CUDA_SDK?=/home/user/NVIDIA_GPU_Computing_SDK3.2
CUDA_TK?=/usr/local/cuda
CUDA_INC=-I$(CUDA_SDK)/C/common/inc -I$(CUDA_TK)/include
CUDA_LIBS=-L$(CUDA_TK)/lib64 -lcudart

GPUCC = $(CUDA_TK)/bin/nvcc
CC=mpicxx

# The location (and name) of the BLAS/Lapack libraries
BLAS_PATH ?= /usr
BLAS_INCLUDE ?= -I$(BLAS_PATH)/include
BLAS_LAPACK_LIB ?= -L$(BLAS_PATH)/lib -llapack -lf77blas -lcblas -latlas 

# Location of GSL header files and libraries
GSL_PATH ?= /usr/
GSL_INCLUDE ?= -I$(GSL_PATH)/include
GSL_LIB = -L$(GSL_PATH)/lib -lgsl

LIBS=$(BLAS_LAPACK_LIB) $(GSL_LIB) -lcublas $(CUDA_LIBS)

INCLUDES = $(CUDA_INC) $(GSL_INCLUDE) $(BLAS_INCLUDE)

CPU_SRC = reference_glm.cpp tvUtil.cpp fortran_matrix.cpp glm.cpp print_matrix.cpp svd.cpp
GPU_SRC = plm.cu
HEADERS= type.h fortran_matrix.h print_matrix.h glm.h plm.h fortran_matrix.h print_matrix.h svd.h cuda_blas.cu
SRC=$(CPU_SRC) $(GPU_SRC)
target = reference_glm
objects = $(patsubst %.cpp,%.o,$(CPU_SRC)) $(patsubst %.cu,%.o,$(GPU_SRC))
all: $(target) convertToBinary

convertToBinary: convertToBinary.o

clean:
	rm -f *~ *.P *.o $(target) convertToBinary

$(target): $(objects)
	$(CC) -o $@ $^ $(LIBS)

%.o:%.cpp
	$(CC) $(DBG) $(OPT_FLAGS) $(INCLUDES) $^ -c
%.o:%.cu
	$(GPUCC) $(DBG) $(OPT_FLAGS) $(INCLUDES) $(CUDA_FLAGS) $^ -c
