#VARIABLES = definition 
F90 := gfortran

SRC_FILES := $(wildcard *.f90)
OBJ_FILES := $(patsubst %.f90, %.o, $(SRC_FILES))

# definitions for compiling
all: $(SRC_FILES)
	$(F90) -fopenmp -Ofast -g -c $^
	$(F90) -fopenmp -g -o PolyMix $(OBJ_FILES)

clean:
	rm -f *.o PolyMix
