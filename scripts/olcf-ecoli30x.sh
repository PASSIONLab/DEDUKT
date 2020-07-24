#!/bin/bash

DATA=ecoli30x
DATA_FILES=all_fastq.txt
DATA_DIR=$PROJWORK/csc103/mme/pacbio_${DATA}
K=17
H=7
CMP_MAX=8

link_files() {
  ln -s ${DATA_DIR}/${DATA_FILES}
  while read NAME; do
    ln -s ${DATA_DIR}/${NAME}
  done < ${DATA_FILES}
}