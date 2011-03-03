# this project seems to require g++ 4.1
CC = g++-4.1

# The location (and name) of the BLAS/Lapack libraries
BLAS_LAPACK_LIB = -lf77blas -latlas -L/usr/lib/atlas/ -llapack

# Location of GSL header files and libraries
#GSL_INCLUDE = /sw/include
GSL_LIB = -lgsl

LIBS=$(BLAS_LAPACK_LIB) $(GSL_LIB)
# You should not need to edit anything below this line...

# Optimization flags (switch this to -O0 -g for debugging)
OPT_FLAGS = -g

# Set all include flags here for later use
#INCLUDE_FLAGS = $(GSL_INCLUDE)

# List of all files in this directory which I can execute
#EXEC_FILES = $(shell find . -type f -perm -u+x -maxdepth 1)

# Command to produce dependencies (.d files) from sources
MAKEDEPEND = $(CC) -M -o $*.d $< -c


SRCS = reference_glm.C
targets = $(patsubst %.C,%,$(SRCS))

all: $(targets)

clean:
	rm -f $(EXEC_FILES) *~ *.P *.o $(targets)

# Override built-in rule for %.o from %.C.  Define this
# and leave it blank!  Without this, Make will try to compile
# directly from %.C to executable without using our stuff below...
% : %.C


# A make rule for automatically generating dependencies,
# http://mad-scientist.net/make/autodep.html
%.o : %.C
	$(MAKEDEPEND); \
	cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
            -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d
	$(CC) $(OPT_FLAGS) $< -o $@ -c


# Note: Running 'make rank' for example, will delete the .o file as an
# itermediate.  You can avoid this behavior by making the .o files "PRECIOUS":

# .PRECIOUS: $(SRCS:.C=.o)

% : %.o
	$(CC) -o $@ $< $(LIBS)

-include $(SRCS:.C=.P)

