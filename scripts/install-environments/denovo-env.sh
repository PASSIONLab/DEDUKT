#!/bin/bash

module load PrgEnv-gnu/7.1
module load cmake/3.7.2 # need 3.5 or higher to support set(CMAKE_CXX_STANDARD 11) with icpc
module load openmpi/1.8.7

CC=$(which mpicc)
CXX=$(which mpic++)

INSTALL_DIR=$BSCRATCH/diBELLA-install
