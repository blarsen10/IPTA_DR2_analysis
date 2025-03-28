#!/bin/bash
#SBATCH --job-name=lite_uf43_CRN13_g4p3_advnoise_20250215
#SBATCH --partition=scavenge
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/lite_uf43_CRN13_g4p3_advnoise_20250215.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/lite_uf43_CRN13_g4p3_advnoise_20250215.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --requeue
#SBATCH --time=1-
#SBATCH --array=0-8
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

Niter=800000
gwb_Nfreqs=13
dataset=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/lite_unfiltered_43_ePSRs
outdir=/vast/palmer/scratch/mingarelli/bbl29/IPTA_DR2_analysis/chains/dr2lite_unfiltered_43/CRN${gwb_Nfreqs}_g4p3_advnoise
emp_dist_path=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/empdists/dr2lite_unfiltered_advnoise.pkl
noisedict_path=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/noisedicts/dr2lite_unfiltered_advnoise.json


srun -n $SLURM_NTASKS python3 ./scripts/run_gwb_adv_noise.py --fixed_gamma --gwb_Nfreqs $gwb_Nfreqs --dataset $dataset --Niter $Niter --outdir $outdir --emp_dist_path $emp_dist_path --noisedict_path $noisedict_path --human bjorn.larsen

