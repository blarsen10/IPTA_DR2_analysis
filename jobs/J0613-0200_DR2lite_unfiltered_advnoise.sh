#!/bin/bash
#SBATCH --job-name=DR2lite_unfiltered_J0613-0200_advnoise_20250213
#SBATCH --partition=scavenge
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/J0613-0200_DR2lite_unfiltered_advnoise_20250213.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/J0613-0200_DR2lite_unfiltered_advnoise_20250213.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-
#SBATCH --array=0-9
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu
#SBATCH --requeue

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

psrname=J0613-0200
dataset=lite_unfiltered_43
Niter=350000
outdir=/vast/palmer/scratch/mingarelli/bbl29/IPTA_DR2_analysis/chains/dr2lite_unfiltered/advnoise/J0613-0200
psr=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/lite_unfiltered_43_ePSRs/J0613-0200.pkl
emp_dist_path=


srun -n $SLURM_NTASKS python3 ./scripts/run_custom_noise.py --psrname $psrname --Niter $Niter --outdir $outdir --psr $psr --emp_dist_path $emp_dist_path --human bjorn.larsen --dataset $dataset

