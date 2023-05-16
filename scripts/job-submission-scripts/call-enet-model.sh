#!/bin/bash
#COBALT -n 1
#COBALT -t 01:00:00
#COBALT -A covid-ct
#COBALT -q debug-flat-quad
#COBALT -o /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/enet_model_log.out
#COBALT -e /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/enet_model_log.err
#COBALT --attrs mcdram=cache:numa=quad:filesystems=home,theta-fs0

echo "COBALT_JOBID = " $COBALT_JOBID
echo "COBALT_JOBSIZE (nodes) =" $COBALT_JOBSIZE
echo "COBALT_PARTNAME = " $COBALT_PARTNAME

NNODES=${COBALT_JOBSIZE}
NRANKS_PER_NODE=1
NTHREADS_PER_CORE=2
NDEPTH=64

NTOTRANKS=1 #$(( NNODES * NRANKS_PER_NODE ))

aprun -n ${NTOTRANKS} -N ${NRANKS_PER_NODE} -d ${NDEPTH} -j ${NTHREADS_PER_CORE} -cc depth /projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/call-enet-model.sh

status=$?

echo "Exit status of aprun is: $status"
exit $status