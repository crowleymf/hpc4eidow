# Model Folder

This folder holds the simulation configurations of different runs


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
