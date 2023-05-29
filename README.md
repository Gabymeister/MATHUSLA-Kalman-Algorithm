# MATHUSLA-Tracker
Tracking software for MATHUSLA Experiment

## Introduction
This repository, including this documentation, is a modification of the following code from: https://github.com/seg188/MATHUSLA-MLTracker 

This repository containts a Kalman Filter Tracking Algorithm based on architecture built for the Linear Tracker above.

## Installation

Most of the dependencies and necessary libraries for this code are available through the anaconda software. If you do not have anaconda, it can be installed following the directions here: https://docs.anaconda.com/anaconda/install/

Once anaconda is installed, the enviornment can be created using the .yml file available in this repository as follows. From the top directory:

```bash
$ conda env create -f env/environment.yml
$ conda activate tracker
```
You will also need the Eigen Library for c++ and to install the joblib package to python.

Note: joblib can be installed through anaconda, Eigen does not require installation. You can downlowd the tar ball from here https://eigen.tuxfamily.org/dox/index.html , then unzip it. If cmake cannot find the directory automatically, you will need to point it to the unzipped eigen directory by setting the EIGEN3_INCLUDE_DIR variable in /tracker/CMakeLists.txt and commenting out the corresponding include_directories command. 

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

Tracker parameters are stored in a txt file located at `tracker/run/par_card.txt`
The parameters are explained in the table below:

| parameter name | usage | unit|
|:--------------|:-------------------------|:---|
|branch              | 0 or 1, 0 for nomal mode and 1 for COSMIC mode |
|debug               | 0 of 1, 1 to turn on debug information| |
|seed                | (float, default =1.0) initial seed value for random generator. Set to -1 to use arbiturary seed| |
|seed_interval               | (float, default =1.0) Maximum interval for track seeding. Interval defined as ds^2 = dr^2-(c*dt)^2| |
|kalman_chi_s                | DEPRECATED                   | |
|kalman_chi_add              | (float, default=200) The maximum accepted chi2 increment for new hit added to the Kalman filter                   | |
|kalman_track_chi                | (float, default=15) Cut on final track reduced-chi2 after smoothing| |
|kalman_pval_drop                | (float, default=1.0) Cut on the smoothed chi2 during dropping steps. If the P value when a hit is added is larger than this number, the hit is dropped                 | |
|kalman_pval_add             | 0.99                 | |
|kalman_pval_track               | 0.95                 | |
|p               | 500.0                    | |
|scint_efficiency                | 0.001                    | |
|merge_cos_theta             | -2.0                 | |
|merge_distance              | 0.0                  | |
|seed_closest_approach               | 300.0                    | |
|vertex_chi2             | 15.0                 | |
|closest_approach_add                | 150.0                    | |
|kalman_vertex_chi_add               | 100000.0                 | |
|kalman_vertex_chi               | 100.0                    | |
|kalman_v_add[0]                | 0.8                   | |
|kalman_v_add[1]                | 1.2                   | |
|kalman_v_drop[0]               | -9999                 | |
|kalman_v_drop[1]               | 9999                  | |
|start_ev                | 0.0                  | |
|end_ev              | 200000.0                 | |
|noise_hz                | 0.0                  | |



