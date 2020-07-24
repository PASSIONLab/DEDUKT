#!/bin/bash
module swap PrgEnv-pgi PrgEnv-gnu
module swap gcc gcc/7.1.0
module load cmake3/3.6.1

CC=$(which cc)
CXX=$(which CC)

PROJECT=bif115 #csc103 # TODO: generalize / add INSTALL_DIR override
INSTALL_DIR=$MEMBERWORK/$PROJECT/diBELLA-install-$(git rev-parse --short HEAD)

