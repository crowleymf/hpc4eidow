# Model Folder

This folder holds the simulation configurations of different runs

# Adding a run
In order to add a run copy the "debug" folder in the same directory with a new unique name describing the simulations

Then update the config.in and data.in to change the simulation parameters; you can also replace model.out.zip with a zip of the file you would like to use (please keep the same name).

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
