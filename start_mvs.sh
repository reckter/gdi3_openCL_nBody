#!/bin/sh
#BSUB -J GDI3_P5_GruppenNR<XY>
#BSUB -W 00:05
#BSUB -x
#BSUB -o /home/kurs2/<TU-ID>/GDI3/Test.%J.out
#BSUB -n 1
#BSUB -q kurs2 

cd /home/kurs2/<TU-ID>/GDI3/

mkdir -p $LSB_BATCH_JID

cd $LSB_BATCH_JID

env | tee -a myLogfile.log

hostname | tee -a myLogfile.log

module load cuda
module load gcc

cd /home/kurs2/<TU-ID>/GDI3/
/home/kurs2/<TU-ID>/GDI3/nBody -k 1 | tee -a myTest -i
