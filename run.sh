#!/bin/bash
#SBATCH --partition=standard --output=output_%a.txt -c 1 -t 5-00:00:00 --mem-per-cpu=190gb
#SBATCH -a 0-29
python lasspia.py configs/cmassS_1_coarse.py routines/combinatorial.py --iJob $SLURM_ARRAY_TASK_ID --nJobs 30

