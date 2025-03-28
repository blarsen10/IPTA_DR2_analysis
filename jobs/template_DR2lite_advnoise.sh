#!/bin/bash
#SBATCH --job-name=DR2lite_{{PSRNAME}}_advnoise_{{DATE}}
#SBATCH --partition={{PARTITION}}
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/{{PSRNAME}}_DR2lite_advnoise_{{DATE}}.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/{{PSRNAME}}_DR2lite_advnoise_{{DATE}}.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time={{RUNTIME_DAYS}}-
#SBATCH --array=0-{{ARRAY_MAX}}
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu
#SBATCH --requeue

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

psrname={{PSRNAME}}
dataset=lite
Niter=350000
outdir=/vast/palmer/scratch/mingarelli/bbl29/IPTA_DR2_analysis/chains/dr2lite/advnoise/{{PSRNAME}}
psr={{PSR}}
emp_dist_path=


srun -n $SLURM_NTASKS python3 ./scripts/run_custom_noise.py --psrname $psrname --Niter $Niter --outdir $outdir --psr $psr --emp_dist_path $emp_dist_path --human bjorn.larsen --dataset $dataset

