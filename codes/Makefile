#VARIABLES = definition 
F90 := gfortran

SRC_FILES := $(wildcard *.f90)
OBJ_FILES := $(patsubst %.f90, %.o, $(SRC_FILES))
MOD_FILES := param.mod correlation.mod pmtypes_module.o pm_subs.o
OPT := -Ofast

# definitions for compiling
all: PolyMix

PolyMix: $(OBJ_FILES) 
	$(F90) -fopenmp -o PolyMix $(OBJ_FILES)

PolyMix.o: param.mod correlation.mod $(MOD_FILES) chaindynamics.o chaincalcs.o PolyMix.f90
	$(F90) -fopenmp $(OPT) -c PolyMix.f90

autocorrelate.o: autocorrelate.f90 param.mod
	$(F90) -fopenmp $(OPT) -c autocorrelate.f90

chaindynamics.o: chaindynamics.f90 correlation.mod 
	$(F90) -fopenmp $(OPT) -c $<

biasd.o: biasd.f90 pm_subs.o
	$(F90) -fopenmp $(OPT) -c biasd.f90

boxcalcs.o: boxcalcs.f90
	$(F90) -fopenmp $(OPT) -c $^

chaincalcs.o: chaincalcs.f90
	$(F90) -fopenmp $(OPT) -c $^

vel.o: vel.f90
	$(F90) -fopenmp $(OPT) -c $^

pm_subs.o: pm_subs.f90 pmtypes_module.o chaincalcs.o chaindynamics.o
	$(F90) -fopenmp $(OPT) -c pm_subs.f90


correlation.mod: corr_module.o
corr_module.o:   corr_module.f90 param.mod
	$(F90) -fopenmp $(OPT) -c corr_module.f90

param.mod: param_module.o
param_module.o: param_module.f90 pmtypes.mod
	$(F90) -fopenmp $(OPT) -c param_module.f90

pmtypes.mod:      pmtypes_module.o
pmtypes_module.o: pmtypes_module.f90
	$(F90) -fopenmp $(OPT) -c $^

clean:
	rm -f *.o PolyMix *.mod
