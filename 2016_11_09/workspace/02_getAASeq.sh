#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -p defq
#SBATCH --array=0-167
#SBATCH -t 5-00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBTACH --mem=8000
cd /home/vxue/data/experimental/SORTCERY/2016_11_09/workspace

if [ ! -d seqframe ]; then
    mkdir seqframe;
fi

if [ ! -d dnaframe ]; then
    mkdir dnaframe;
fi

if [ ! -d qualframe ]; then
    mkdir qualframe;
fi


source activate localEnv2
python getDNASeq_mp.py ${SLURM_ARRAY_TASK_ID}

