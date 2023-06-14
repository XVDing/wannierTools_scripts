#!/bin/bash
#BSUB -q sr850   #  sr850 sr860-768 sr860-1536
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -J dxy
#BSUB -R "span[ptile=40]"
#BSUB -R "select[hname!='r13n18']"
hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
#-------------intelmpi+ifort------------------------------------------
source /work/software/intel/bin/compilervars.sh -arch intel64 -platform linux
#---------------------------------------------------------------------

#~/dxy/software/miniconda/py39/bin/python findnodes.py > nodes.log

#~/dxy/software/miniconda/py39/bin/python -u nodes_arc.py > arc_nodes.log

~/dxy/software/miniconda/py39/bin/python -u /work/wangr/data/dxy/scripts/program/wt/transform_nodes-firstBZ.py > trans.log
