#!/bin/bash
module swap PrgEnv-intel PrgEnv-gnu # TODO make compatible for any PE/compiler set -- must satisfy (__cplusplus < 201300) - may want to change code seqan/platform.h ln 154
module unload cmake #used instead of swap so that whether or not cmake is loaded the next line will succeed
module load  cmake/3.14.4

CC=$(which cc)
CXX=$(which CC)

INSTALL_DIR=$SCRATCH/diBELLA-install-$(git rev-parse --short HEAD)

if [ ! -z "$KNL" ]; then
  echo "Building for Cori KNL nodes"
  INSTALL_DIR=$SCRATCH/diBELLA-install-knl-$(git rev-parse --short HEAD)
  module swap craype-haswell craype-mic-knl
  DEFS+=" -DCFLAGS=-march=knl -DCXXFLAGS='-march=knl -ldmapp'"
  export MPICH_USE_DMAPP_COLL=1
  export MPICH_NETWORK_BUFFER_COLL=1
fi
