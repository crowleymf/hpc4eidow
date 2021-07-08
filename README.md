[![Build, run and test](https://github.com/jamilgafur/hpc4eidow/actions/workflows/fortran.yml/badge.svg)](https://github.com/jamilgafur/hpc4eidow/actions/workflows/fortran.yml)


# Allocating an Interactive Node on Eagle

```bash
srun --time=1:00:00 --account=hpc4eidow -ntasks=1 --partition=debug --pty $SHELL
```
# Submitting the code to run

```bash
sbatch srun.pbs
```
# Running the code on an Interactive Node
```batch
# allocate an interactive node
./inter_srun

# loads the gfortran compiler
module load gcc

# run the module interactivly
source run.pbd
```

# Directory structure

```text
─── codes: holds the actual simulation codes
│   ├── PolyMix.f90: Main program that implements the Cooperative Motion Algorithm
│   │   developed by Tadeusz Pakula;This program is for bipolar shear flow of linear chains.
│   ├── autocorrelate.f90: ???
│   ├── biasd.f90: makes Monte Carlo random direction moves for each particle
│   ├── boxcalcs.f90: calculates the stress, lattice site density and segmental desntiy
│   ├── build: compiles and builds the executable ``PolyMix''
│   ├── chaincalcs.f90: Calculates the end-to-end vector, radius of gyration, order
│   │     parameters and center of mass density
│   ├── chaindynamic.f90: The chaindynamics subroutine presented here is used to
│   │     calculate the mean squared displacement variables called g1, g2, g3, g4, g5.
│   │     These variables track the motion of of various monomers. From these
│   │     correlation functions we can obtain the dynamic "finger print" of the system
│   └── vel.f90: generates velocity profiles
├── model: holds the input data to build, and run the simulation according to user defined specifications
│   ├── conf.in: the configuration file for the model
│   ├── data.in: the data file for the model
│   ├── inter_srun: code used to allocate an interactive node on eagle
│   ├── readme: temp readme
│   └── run.pbs: runs the simulation
├── README.md: Head Readme File
└── .gitignore: files to ignore

```
