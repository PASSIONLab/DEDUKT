#!/bin/bash
#
QUEUE=debug
ACCT="m2865"
LICENSE="SCRATCH"
TIME="00:30:00"
MIN_TIME="00:00:30"
NAME="dibella-single"
NODES=1
PPN=24
CPP=1
DBG="-dbg" #""
STAGES=""
SCRIPT="run_dibella.slurm"
LINKER="$CSCRATCH/pacbio_ecoli30x/linker.sh" #<-- update
OTHER="--test-only" #"-d afterany:9777467"

sbatch -A ${ACCT} -L ${LICENSE} -J ${NAME} -q ${QUEUE} ${OTHER} -t ${TIME} --time-min ${MIN_TIME} -N ${NODES} --export=PPN=${PPN},CPP=${CPP},DBG=${DBG},STAGES="$STAGES",LINKER=$LINKER ${SCRIPT} && date

