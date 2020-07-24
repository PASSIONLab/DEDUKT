#!/bin/bash

echo "Setting up build environment:"
ENVSCRIPTS_PATH="scripts/install-environments"
if [ ! -z "$1" ]; then 
  echo "using provided environment script: "$1
  source $1
elif [[ $NERSC_HOST = "edison" ]]; then
  source $ENVSCRIPTS_PATH/edison-env.sh
elif [[ $NERSC_HOST = "cori" ]]; then
  source $ENVSCRIPTS_PATH/cori-env.sh
elif [[ $NERSC_HOST = "denovo" ]]; then
  source $ENVSCRIPTS_PATH/denovo-env.sh
elif [[ $HOST = "titan"* ]]; then
  source $ENVSCRIPTS_PATH/titan-env.sh
elif [[ $(uname -s) = "Darwin" ]]; then			# Mac 
  source $ENVSCRIPTS_PATH/mac-env.sh
fi

# specify defaults if build & install directories and compilers are not provided in the environment
if [ -z "$BUILD_DIR" ]; then
  BUILD_DIR="build"
fi
if [ -z "$INSTALL_DIR" ]; then
  INSTALL_DIR="bin"
  # INSTALL_DIR=$HOME/dibella_bsp/build
fi
if [ -z "$CC" ]; then
  CC=$(which gcc)
fi
if [ -z "$CXX" ]; then
  CXX=$(which g++)
fi

DEFS+=" -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DHIPMER_FULL_BUILD=1 "

# setup debug/release
if [ $DEBUG ]; then
  BUILD_DIR=$BUILD_DIR"-dbg"
  INSTALL_DIR=${INSTALL_DIR}"-dbg"
  DEFS+=" -DCMAKE_BUILD_TYPE=Debug "
else
  DEFS+=" -DCMAKE_BUILD_TYPE=Release "
fi

DEFS+=" -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR "
echo $DEFS

mkdir $BUILD_DIR
RETURN_DIR=$(pwd)
cd $BUILD_DIR
VERBOSE=1 cmake $DEFS ..

echo "Install directory is $INSTALL_DIR"
# read -n 1 -p "Press 'y' to install:" yes_install
# if [ "$yes_install" = "y" ]; then
  VERBOSE=1 make all install
# fi

cd $RETURN_DIR
