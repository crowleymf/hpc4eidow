#VARIABLES = definition 
F90 := gfortran

SRC_FILES := $(wildcard *.f90)
OBJ_FILES := $(patsubst %.f90, %.o, $(SRC_FILES))
OPT := -Ofast

# definitions for compiling
all: PolyMix

PolyMix: $(OBJ_FILES) 
	$(F90) -fopenmp -o PolyMix $(OBJ_FILES)

PolyMix.o: PolyMix.f90 pmtypes_module.o param_module.o corr_module.o chaindynamics.o chaincalcs.o boxcalcs.o autocorrelate.o
	$(F90) -fopenmp $(OPT) -c PolyMix.f90

autocorrelate.o: autocorrelate.f90 param_module.o pm_subs.o pmtypes_module.o 
	$(F90) -fopenmp $(OPT) -c autocorrelate.f90

chaindynamics.o: chaindynamics.f90 corr_module.o pmtypes_module.o
	$(F90) -fopenmp $(OPT) -c $<

boxcalcs.o: boxcalcs.f90 pmtypes_module.o param_module.o pm_subs.o
	$(F90) -fopenmp $(OPT) -c boxcalcs.f90

chaincalcs.o: chaincalcs.f90 pmtypes_module.o param_module.o
	$(F90) -fopenmp $(OPT) -c chaincalcs.f90

vel.o: vel.f90 pmtypes_module.f90 param_module.o
	$(F90) -fopenmp $(OPT) -c vel.f90

pm_subs.o: pm_subs.f90 pmtypes_module.o param_module.o
	$(F90) -fopenmp $(OPT) -c pm_subs.f90


corr_module.o:   corr_module.f90 param_module.o pmtypes_module.o
	$(F90) -fopenmp $(OPT) -c corr_module.f90

param_module.o: param_module.f90 pmtypes_module.o
	$(F90) -fopenmp $(OPT) -c param_module.f90 pmtypes_module.f90

pmtypes_module.o: pmtypes_module.f90
	$(F90) -fopenmp $(OPT) -c $^

clean:
	rm -f *.o PolyMix *.mod
