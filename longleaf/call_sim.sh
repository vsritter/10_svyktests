#!/bin/bash

#SBATCH --job-name=sim_svytest
#SBATCH --output=./longleaf/slurm_out/slurm_log_%a.out

#SBATCH --array=1-20
#SBATCH --ntasks=24
#SBATCH --partition=general
#SBATCH --time=48:00:00
#SBATCH --mem=50g

#SBATCH --mail-type=end
#SBATCH --mail-user=vritter@email.unc.edu

# Create output dir for this job
OUTPUT_DIR=./longleaf/slurm_out/job_${SLURM_ARRAY_JOB_ID}
mkdir -p ${OUTPUT_DIR}

# Source code dir
SOURCE_FILE=./longleaf/call_sim

# Record initial time
TIME=$(date +%m-%d-%Y--%H:%M:%S)

# Call R script
R CMD BATCH --no-save --no-restore ${SOURCE_FILE}.R ${OUTPUT_DIR}/r_log_task_${SLURM_ARRAY_TASK_ID}.Rout

# Move slurm logs to job folder (silently)
mv ./longleaf/slurm_out/*.out ${OUTPUT_DIR} 2> /dev/null


if [ ${SLURM_ARRAY_TASK_ID} == 1 ]
then
  echo $TIME '|' $SLURM_ARRAY_JOB_ID '|' $ENV_REP'|' $ENV_TEST'|' $ENV_SIMID'|' $ENV_SEED'|' $ENV_S_SIZE'|' $ENV_CENS'|' $ENV_PSU_SIZE'|' $ENV_EFFECT '|' 0 | tee -a ./longleaf/slurm_out/log.txt
fi

# happy end
exit 0
