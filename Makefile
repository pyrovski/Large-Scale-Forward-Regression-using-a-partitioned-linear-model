ifeq ($(dbg),1)
DBG=-D_DEBUG
OPT_FLAGS = -g
else
DBG=
OPT_FLAGS = -O3
endif


CUDA_SDK=/home/user/NVIDIA_GPU_Computing_SDK3.2
CUDA_TK=/usr/local/cuda
CUDA_INC=-I$(CUDA_SDK)/C/common/inc -I$(CUDA_TK)/include

CC = $(CUDA_TK)/bin/nvcc

# The location (and name) of the BLAS/Lapack libraries
BLAS_LAPACK_LIB = -lf77blas -latlas -L/usr/lib/atlas/ -llapack
# You should not need to edit anything below this line...

# Location of GSL header files and libraries
#GSL_INCLUDE = /sw/include
GSL_LIB = -lgsl

LIBS=$(BLAS_LAPACK_LIB) $(GSL_LIB)

# Set all include flags here for later use
#INCLUDE_FLAGS = $(GSL_INCLUDE)

# List of all files in this directory which I can execute
#EXEC_FILES = $(shell find . -type f -perm -u+x -maxdepth 1)

# Command to produce dependencies (.d files) from sources
MAKEDEPEND = $(CC) $(DBG) -M -o $*.d $< -c $(CUDA_INC)

SRCS = reference_glm.cu
targets = $(patsubst %.cu,%,$(SRCS))

all: $(targets)

clean:
	rm -f $(EXEC_FILES) *~ *.P *.o $(targets)

# Override built-in rule for %.o from %.cpp.  Define this
# and leave it blank!  Without this, Make will try to compile
# directly from %.cpp to executable without using our stuff below...
% : %.cu


%.o : %.cu
	$(CC) $(DBG) $(OPT_FLAGS) $(CUDA_INC) $< -o $@ -c



% : %.o
	$(CC) -o $@ $< $(LIBS)


