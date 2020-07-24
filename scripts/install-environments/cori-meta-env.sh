#!/bin/bash
module swap PrgEnv-intel PrgEnv-gnu # TODO make compatible for any PE/compiler set -- must satisfy (__cplusplus < 201300) - may want to change code seqan/platform.h ln 154
module unload cmake #used instead of swap so that whether or not cmake is loaded the next line will succeed
module load cmake/3.8.2 # need 3.5 or higher to support set(CMAKE_CXX_STANDARD 11) with icpc

CC=$(which cc)
CXX=$(which CC)
#TODO set reasonable values for metagenomes w.r.t. MAX_NUM_READS_LIST
CMAKE_FLAGS=$CMAKE_FLAGS" -DTIGHT_READS=1 -DKMERA_ONLY=1"
BUILD_DIR=build-$USER-$(git rev-parse --short HEAD)
INSTALL_DIR=$SCRATCH/metaDiBELLA-$(git rev-parse --short HEAD)-install

