# Colony image analysis of Oa and a collection of 28 partners in presence and absence of thiamin

Brief description of the project and what the code in this repository is used for.

## 1. Repository contents

This repository contains the code used to process and analyze the data associated with this project. The raw and processed data are stored separately and are **not included in this repository**.

The repository is organized as follows:

```text
repository/
├── analysis/
│   ├── ...
│   └── ...
├── scripts/
│   ├── ...
│   └── ...
├── environment.yml
└── README.md
```

## 2. Data availability and organization

The data used by the analysis scripts are stored separately from this Git repository.

> The data is stored on the NAS under nasdcsr/FAC/FBM/DMF/smitri/nccr_microbiomes/D2c/colony_patterns/.

The data stored on the NAS are organized as follows:

```text
data/
├── raw/
│   ├── experiment_1/
│   ├── experiment_2/
│   └── ...
├── processed/
│   ├── ...
│   └── ...
└── ...
```

### Raw data

`data/raw/` contains the original experimental data.

[Describe the important file types and naming conventions.]

For example:

```text
experiment_1/
├── image_001.tif
├── image_002.tif
└── ...
```

### Processed data

`data/processed/` contains files generated from the raw data and used as input for subsequent analyses.

[Describe what processing was performed and which scripts generate these files.]

For example, microscopy images were classified using ilastik 1.4.0b21. The resulting probability maps (`.h5`) are used as input for the image-analysis scripts.

## Connecting the code to the data

The data directory is not part of this repository. Before running the analysis, the user must specify the location of the data on their own system.

[CHOOSE THE APPROACH ACTUALLY USED BY THE CODE.]

For example, the project may have a configuration file containing:

```text
DATA_DIR=/path/to/data
```

or a Python variable:

```python
DATA_DIR = "/path/to/data"
```

Replace this path with the location of the data directory on your system.

The analysis scripts then construct paths relative to `DATA_DIR`, rather than relying on absolute paths specific to the original computer.

For example:

```text
DATA_DIR
└── raw/
    └── experiment_1/
```

corresponds to:

```python
DATA_DIR / "raw" / "experiment_1"
```

## Analysis workflow

The overall workflow is:

```text
Raw experimental data
        │
        ▼
Preprocessing / ilastik classification
        │
        ▼
Processed data / probability maps
        │
        ▼
Image-analysis scripts
        │
        ▼
Quantified measurements
        │
        ▼
Statistical analysis and figures
```

[Optionally list the scripts responsible for each step.]

| Step | Input | Script/software | Output |
|---|---|---|---|
| 1. Pixel classification | Raw microscopy images | ilastik 1.4.0b21 | `.h5` probability maps |
| 2. Image analysis | Probability maps | `script_name.py` | Quantified colony properties |
| 3. Statistical analysis | Quantified data | `script_name.py` | Statistics and figures |

## Software environment

The Python environment required for the analysis can be created using:

```bash
conda env create -f environment.yml
conda activate colony_patterns
```

Pixel classification was performed separately using **ilastik 1.4.0b21**.

## Reproducing the analysis

To reproduce the analysis:

1. Clone this repository.
2. Download or obtain the associated data.
3. Organize the data according to the directory structure described above.
4. Set the data path to the appropriate location on your system.
5. Create and activate the Conda environment.
6. Run the analysis scripts in the order described under **Analysis workflow**.

The analysis code uses relative paths derived from the configured data directory so that the location of the data can differ between systems.
Scripts developed by Alessia Del Panta at the University of Lausanne, Department of Fundamental Microbiology.

Jupyter notebook tested on Ubuntu 22.04 on an Intel Core i10 processor under Python 3.7.11 using jupyter-core 4.9.2. !! CHECK !!

### Pixel classification

Pixel classification was performed using ilastik 1.4.0b21 with the Pixel
Classification workflow. The trained ilastik project (`.ilp`) used for
classification is provided in [path].

Probability maps exported from ilastik (`.h5`) are provided in [path] and
were used as input for the subsequent image-analysis pipeline.

### Instructions for Installation

1. Clone github repository !! INSTRUCTIONS !!
2. Assuming conda is installed, install the environments and dependencies using 
`conda env create -f colony_patterns.yml` and
`conda env create -f colony_patterns_glm.yml`
Canghe the path ????
4. Activate the environment by typing `conda activate mir430-distance-analysis`
5. Launch a jupyter instance by typing `jupyter notebook`

### Instructions for use 

1. In a jupyter notebook, open the `miR430-distance-analysis.ipynb` notebook, and execute. 
2. The output of this script is provided in the form of four files:
    - `miR430-distance-dataframe.csv` contains the trajectories of the distances between the tracked loci.
    - `miR430-speed-time-dataframe.csv` contains, for each trajectory, the duration of each of the condensed states, and the time it takes to go from the de-compacted to the condensed state, for each of the compacted states. 
    - `miR430-summary-stats.csv` contains summary statistics such as the number of oscillations, or the average time in condensed state, separated into Nanog present/absent conditions. 
    - `hp_tests.csv` contains the p-values for the Mann-Whitney tests performed to understand the role of Nanog in the dynamics of the loci. 

### Running the software on your data

This software was written to operate on the data generated for the current publication. In order to run this analysis on your own data, you must provide a file with the trajectories of miR430 loci. The header should be the following header:

Date	Stage	Time	TrackID	Position.X	Position.Y	Position.Z	Nucleus.number	Condition

Where:

- "Date" is the date at which the loci where imaged
- "Stage" is the stage at which the cell was at during imaging
- "TrackID" is a unique ID for each allele within a nucleus
- "Time" is an integer indicating the time step at which a locus was imaged
- "Position.X", "Position.Y", "Position.Z" are the coordinates of an imaged lcus at a certain time step
- "Nucleus.number" is a unique ID for each nucleus
- "Condition" is either "Mutant" or "Inhibition"

Then, run the script, updating the file path (variable "filename") and the time interval between frames (variable "delta_t"). 
