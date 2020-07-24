#!/bin/bash
#
# Variables (expected from $1)		Description
# 	link_files()                    function for linking data files to run directory
# 	K                               k-mer length
# 	H                               reliable k-mer upper limit
#
# Variables (set externally) 		Options
#	DBG                             "-dbg", ""
#	INMEMORY                        "-inmemory", ""
#	PPN	(processors per node)	integer [1,16] on Titan)
#	TASKS_PER_NUMA			integer [1, 8] on Titan, must be < PPN
#	CPP (cores per processor)	integer [1, 2] on Titan
#	CMP_MAX (precompiled readId buffer lengths, usually next greater power of 2 after $H)
#
# Variables (internal - may want to change)
#	BIN (name of binary to execute)
#

module swap PrgEnv-intel PrgEnv-gnu
module list
######### Set data-dependent variables #########
if [ -z $LINKER ]; then
	echo "Script for linking data files is required. Set LINKER."
	exit 1
fi
source $LINKER
E=2

######### Set performance/parallelism variables #########
[ -z $PPN ] && PPN=${SLURM_CPUS_ON_NODE}
[ -z $CPP ] && CPP=1
if [ -z LOWER ]; then LOWER=$SLURM_JOB_NUM_NODES; fi
NUM_NODES=$SLURM_JOB_NUM_NODES
while [ $NUM_NODES -ge $LOWER ]; do # begin scaling loop
MPI_TASKS=$(($NUM_NODES * $PPN))
RUN="srun -n ${MPI_TASKS} -c ${CPP} "

######### Set runtime variables #########
INSTALL_DIR=${SCRATCH}/diBELLA-install${DBG}/bin
RUN_DIR=${SCRATCH}/diBELLA${DBG}_${SLURM_JOB_NUM_NODES}n_${PPN}c_${DATA}_k${K}_e${E}_h${H}_${SLURM_JOB_ID}
mkdir ${RUN_DIR}
cd ${RUN_DIR}
link_files                        # function included/required in $1, includes $DATA_FILES, $H, etc.
BIN="main-${CMP_MAX}-32"		
CMD=$RUN" ${INSTALL_DIR}/${BIN} -N $PPN -k ${K} -i ${DATA_FILES} -e ${E} -u ${H} ${STAGES}"

######### Report runtime setup #########
echo "cleaning-out previous log files"
set -x
rm per_thread/*/*/*.log
set +x

######### Report runtime setup #########
BAR="===================="
echo "Host Name: "$HOSTNAME
echo "Number of nodes: $NUM_NODES of $SLURM_JOB_NUM_NODES"
echo "Data: "$DATA
if [ -z ${DBG} ]; then
  echo "Running in mode: Release"
else
  echo "Running in mode: Debug"
fi
echo "Run directory is: "$(pwd)
echo "Date: $(date)"
echo $BAR" Starting diBELLA pipeline "$BAR
echo $CMD

######### Run #########
if [ -z $SCRIPTTEST ]; then $CMD ; fi #> $OUTFILE 2>&1 ; fi # if not a script test, run the code
NUM_NODES=$(($NUM_NODES/2))
done

