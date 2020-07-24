#!/bin/bash
#
QUEUE=regular # debug, regular, or killable
DBG="" 	#"-dbg" or ""
ACCT="m2865"
LICENSE="SCRATCH"
PPN=16
CPN=2
NAME="dibella-scaling"
SCRIPT="run_dibella.slurm"
STAGES="" # all by default
LINKER="" # <-- set
TIME="00:10:00"
MIN_TIME="00:05:00"
NODES=64
OTHER="" #"--test-only" #"-d afterany:9777467"

for i in 1 2 4 8 16 32 64; do
    let N=$NODES/$i
    sbatch -A ${ACCT} -L ${LICENSES} -J ${NAME} -q ${QUEUE} -t ${TIME} --time-min=${MIN_TIME} -N ${N} ${OTHER} --export=PPN=${PPN},CPP=${CPP},DBG=${DBG},STAGES="$STAGES",LINKER=$LINKER ${SCRIPT} && date
done

