#!/bin/bash
module swap PrgEnv-intel PrgEnv-gnu # TODO make compatible for any PE/compiler set -- must satisfy (__cplusplus < 201300) - may want to change code seqan/platform.h ln 154
module unload cmake #used instead of swap so that whether or not cmake is loaded the next line will succeed
module load cmake/3.8.1 # need 3.5 or higher to support set(CMAKE_CXX_STANDARD 11) with icpc

CC=$(which cc)
CXX=$(which CC)

INSTALL_DIR=$SCRATCH/diBELLA-install
