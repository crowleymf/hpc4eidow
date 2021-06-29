# Model Folder

This folder holds the simulation configurations of different runs

# Adding a run
In order to add a run copy the "debug" folder in the same directory with a new
unique name describing the simulations

Then update the config.in and data.in to change the simulation parameters; you
can also replace model.out.zip with a zip of the file you would like to use
(please keep the same name).

## Running your own simulation

Initially you can copy the default run files into a new folder in the same
directory; the format for each folder is as follows

### Conf.in (configurations)
This file contains the FCC lattice configurations in the form of a lookup
table. The first 13 lines represents the current polymers ability to move
to a new location in 3-D space, the inverse of the direction and the
polymer id (between 1 and 13 inclusive). This must be defined for each polymer
position in the configuration.

  x , y , z, -x ,-y,-z , polymer_id

The next set of lines are still being figured output

### data.in

This file defines the input files they are as follows:
```text
1| input model file name
2| output file name
3| (The number of Monte Carlo iteration) (the iteration to write output) (the iteration for logging)
4| -----
5| (the iteration to change the p-value) (the new pvalue)
6| ------
```

# Running multiple simulations in parallel on a single node
within the /mode/run directory you can run the command ```bash make all``` to run all simulations in parallel


# Submitting multiple simulations as batch jobs
within the /mode/run directory you can run the command ```bash make batch``` to submit all simulations to SLURM


# Directory
```text
|── golden: Output files for unit testing
├── runs: Folder to put individual simulation configurations (single simulation per  folder)
│   ├── default: Folder that
|   └── <USER FOLDER>: Users add folders here containing their simulation parameters
├── README.md: This File
└── .gitignore: files to ignore
└── golden_test.py: python unit test to compare the output of runs/default to golden_test
└── inter_srun: allocates an interactive node
```
