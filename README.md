# README

## Analysis script for the dynamics of mir430 loci

Script developed by Alessia Del Panta at the University of Lausanne, Department of Fundamental Microbiology.

Jupyter notebook tested on Ubuntu 22.04 on an Intel Core i10 processor under Python 3.7.11 using jupyter-core 4.9.2.

### Instructions for Installation (typical install time, < 2 mins)

1. Place the `miR430-distance-analysis.ipynb`, `20230613-miR430-spots-raw-data.csv` and `mir430-distance-analysis.yml` files into a folder.
2. Activate a terminal with conda installed and install the environment and dependencies using 
`conda env create -f mir430-distance-analysis.yml`
3. Activate the environment by typing `conda activate mir430-distance-analysis`
4. Launch a jupyter instance by typing `jupyter notebook`

### Instructions for use (typical runtime < 30s)

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
