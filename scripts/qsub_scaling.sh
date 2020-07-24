#!/bin/bash
#
# References:
# [1] https://www.olcf.ornl.gov/for-users/system-user-guides/titan/running-jobs/#titan-batch-queues
# [2] https://www.olcf.ornl.gov/for-users/system-user-guides/titan/running-jobs/#job-priority-by-processor-count
#
# CAN'T REQUEST MORE THAN 12 HRS FOR JOBS WITH <3750 NODES [2] 
# CAN'T REQUEST MORE THAN 6 HRS FOR JOBS WITH <313 NODES [2]
# CAN'T REQUEST MORE THAN 2 HRS FOR JOBS WITH <126 NODES [2]
#
QUEUE=batch 	# debug, batch, or killable
DBG="" 			#"-dbg" or ""
ACCT="CSC103"
INMEMORY=""		#"-inmemory" or ""
PPN=16
CPN=2
NAME="dibella-scaling"
STAGES=""
SCRIPT="run_dibella.pbs"
DATA_SCRIPT=""	# <-- set

WALLTIME="02:00:00"
NODES=64
for i in 1 2 4 8 16 32 64; do
    let N=$NODES/$i
    qsub -S /bin/bash -A ${ACCT} -N ${NAME} -q ${QUEUE} -l walltime=${WALLTIME},nodes=${N} -v PPN=${PPN},CPN=${CPN},DBG=${DBG},STAGES="$STAGES" ${SCRIPT} ${DATA_SCRIPT}
done

