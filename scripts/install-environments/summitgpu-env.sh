#!/bin/bash
module swap xl gcc #/8.1.1
module load cuda
module load cmake/3.15.2

CC=mpicc #gcc #$(which cc)
CXX=mpic++ #g++ #$(which CC)

PROJECT=bif115 #csc103 # TODO: generalize / add INSTALL_DIR override
INSTALL_DIR=$MEMBERWORK/$PROJECT/diBELLA-install-$(git rev-parse --short HEAD)


