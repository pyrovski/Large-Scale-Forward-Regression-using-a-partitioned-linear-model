# 
# 
# Copyright (c) 2011, The Arizona Board of Regents on behalf of 
# The University of Arizona
# 
# All rights reserved.
# 
# Developed by Peter Bailey, Tapasya Patki, and Greg Striemer with
# support from the iPlant Collaborative as a collaboration between
# participants at BIO5 at The University of Arizona (the primary hosting
# institution), Cold Spring Harbor Laboratory, The University of Texas
# at Austin, and individual contributors. Find out more at
# http://www.iplantcollaborative.org/.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are
# met:
# 
#  * Redistributions of source code must retain the above copyright 
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright 
#    notice, this list of conditions and the following disclaimer in the 
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of the iPlant Collaborative, BIO5, The University 
#    of Arizona, Cold Spring Harbor Laboratory, The University of Texas at 
#    Austin, nor the names of other contributors may be used to endorse or 
#    promote products derived from this software without specific prior 
#    written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
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
COMMON_FLAGS+=-D$(CUDA_SM_VER)

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
BLAS_LAPACK_LIB ?= -L$(BLAS_PATH)/lib -llapack -lf77blas -lcblas -latlas -lgfortran

# Location of GSL header files and libraries
GSL_PATH ?= /usr
GSL_LIB_PATH ?= $(GSL_PATH)/lib
GSL_INCLUDE ?= -I$(GSL_PATH)/include
GSL_LIB = -L$(GSL_LIB_PATH) -lgsl

LIBS=$(GSL_LIB) $(BLAS_LAPACK_LIB)

INCLUDES = $(CUDA_INC) $(GSL_INCLUDE) $(BLAS_INCLUDE)

CPU_SRC = reference_plm.cpp tvUtil.cpp fortran_matrix.cpp glm.cpp print_matrix.cpp svd.cpp md5.cpp pinv.cpp
GPU_SRC = plm.cu
HEADERS= type.h fortran_matrix.h print_matrix.h glm.h plm.h fortran_matrix.h print_matrix.h svd.h cuda_blas.cu
SRC=$(CPU_SRC) $(GPU_SRC)
target = reference_plm
objects = $(patsubst %.cpp,%.o,$(CPU_SRC)) $(patsubst %.cu,%.o,$(GPU_SRC))
all: $(target) convertToBinary transpose nonneg

convertToBinary: convertToBinary.o
	$(CC) $(DBG) $(OPT_FLAGS) $(INCLUDES) $(COMMON_FLAGS) $^ -o $@

transpose: transpose.cpp
	$(CC) $(DBG) $(OPT_FLAGS) $(INCLUDES) $(COMMON_FLAGS) $^ -o $@

nonneg: nonneg.cpp
	$(CC) $(DBG) $(OPT_FLAGS) $(INCLUDES) $(COMMON_FLAGS) $^ -o $@

clean:
	rm -f *~ *.P *.o $(target) convertToBinary transpose nonneg

$(target): $(objects)
	$(CC) -o $@ $^ $(LIBS)

%.o:%.cpp
	$(CC) $(DBG) $(OPT_FLAGS) $(INCLUDES) $(COMMON_FLAGS) $^ -c
%.o:%.cu
	$(GPUCC) $(DBG) $(OPT_FLAGS) $(INCLUDES) $(COMMON_FLAGS) $(CUDA_FLAGS) $^ -c
