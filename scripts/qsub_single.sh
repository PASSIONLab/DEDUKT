#!/bin/bash
#
# reference: https://www.olcf.ornl.gov/for-users/system-user-guides/titan/running-jobs/#titan-batch-queues
#
QUEUE=debug # alternative is 'batch'
ACCT="BIF115"
MYPROJECT="bif115"
TIME="00:15:00"
NAME="ecoli30x-dibella"
NODES=32
PPN=16
CPP=2
DBG="" #"-dbg"
STAGES="\"-s all -q 1 -B -a\" " 
SCRIPT="./run_dibella.pbs"
SCALING_LIMIT=1
LINKER="$PROJWORK/$MYPROJECT/pacbio_data/ecoli30x_linker.sh"
INSTALL="/lustre/atlas/scratch/mme/bif115/diBELLA-install-6915300/bin"
echo OLCF
echo $HOST
echo $QUEUE
CMD="qsub -S /bin/bash -A ${ACCT} -N ${NAME} -q ${QUEUE} -l walltime=${TIME},nodes=${NODES} -v PPN=${PPN},CPP=${CPP},DBG=${DBG},LOWER=$SCALING_LIMIT,MYPROJECT=$MYPROJECT,LINKER=$LINKER,INSTALL_DIR=$INSTALL,STAGES="$STAGES" ${SCRIPT}"
#echo $CMD 
$CMD
echo "Queued"
date
echo "dibella"$DBG
(source $LINKER; echo $DATA)
echo $(pwd)

