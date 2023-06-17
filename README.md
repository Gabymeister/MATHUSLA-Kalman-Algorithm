# MATHUSLA-Tracker
Tracking software for MATHUSLA Experiment

## Introduction
This repository, including this documentation, is a modification of the following code from: https://github.com/seg188/MATHUSLA-MLTracker 

This repository containts a Kalman Filter Tracking Algorithm based on architecture built for the Linear Tracker above.

## Installation

The code requires installation of `gcc`,`ROOT`,`eigen` and `dlib` to compile. 

On cedar, `gcc`,`ROOT`,`eigen` can be loaded with the module command, and dlib can be added with an installation in my directory. 
```bash
$ module load qt/5.15.2 gcc/9.3.0 StdEnv/2020   root/6.26.06  eigen/3.3.7
$ PATH=$PATH:/project/def-mdiamond/tomren/mathusla/dlib-19.24/install
```

Now, the project can be built using cmake:

```bash
$ cd tracker
$ cd build
$ cmake ../ 
$ make 
```

At this point, the tracker executable is available in the /build/ directory. Note that the /build/ directory MUST be placed in the /tracker/ directory. 


## Running the Tracker

The tracker requires two command line arguments, the path to an input file, and the path to which the output file should be written. The input file should be the output from a MATHUSLA Mu-Simulation run, a Geant4 based simulation of particle passage through the MATHUSLA detector. 

The Mu-Simulation repository can be found here: https://github.com/MATHUSLA/Mu-Simulation

An example command to run the tracker:

```bash
$ ./tracker path_to_input_file path_to_write_output 
```
A script for automating series runs of the tracker, and further documentation about it, is located in the /run/ directory. 

Job submission scripts for parallelising the tracker or analysis code for large datasets is found in parallel.

### Tracker configuration

Tracker configuration parameters are stored in a txt file located at `tracker/run/par_card.txt`
The parameters are explained in the table below:

| Parameter name | Usage | Default value [unit]|
|:--------------|:-------------------------|:---|
|branch              | 0 or 1, 0 for nomal mode and 1 for COSMIC mode | 0|
|debug               | 0 of 1, 1 to turn on debug information| |
|seed                | initial seed value for random generator. Set to -1 to use arbiturary seed|1.0 |
|seed_interval               |Maximum interval for track seeding. Interval defined as ds^2 = dr^2-(c*dt)^2| 7|
|kalman_chi_s                | DEPRECATED                   | |
|kalman_chi_add              |  The maximum accepted chi2 increment for new hit added to the Kalman filter                   | 200|
|kalman_track_chi                |Cut on final track chi2/ndof after smoothing| 15|
|kalman_pval_drop                | Cut on the smoothed chi2 P-value during dropping steps. If the P value when a hit is added is larger than this number, the hit is dropped                 | 1.0|
|kalman_pval_add             | 0.99                 | |
|kalman_pval_track               | 0.95                 | |
|p               | 500.0                    | |
|scint_efficiency                | 0.001                    | |
|merge_cos_theta             | -2.0                 | |
|merge_distance              | 0.0                  | |
|vertex_seed_dist               | 300.0                    | |
|vertex_chi2             | 15.0                 | |
|vertex_add_max_distance                | 150.0                    | |
|kalman_vertex_chi_add               | 100000.0                 | |
|kalman_vertex_chi               | 100.0                    | |
|kalman_v_add[0]                | 0.8                   | |
|kalman_v_add[1]                | 1.2                   | |
|kalman_v_drop[0]               | -9999                 | |
|kalman_v_drop[1]               | 9999                  | |
|start_ev                | 0.0                  | |
|end_ev              | 200000.0                 | |
|noise_hz                | 0.0                  | |



