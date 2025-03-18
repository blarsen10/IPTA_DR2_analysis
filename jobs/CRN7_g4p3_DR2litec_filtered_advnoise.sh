#!/bin/bash
#SBATCH --job-name=litec_filtered_CRN7_g4p3_advnoise_20250212
#SBATCH --partition=scavenge
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/litec_filtered_CRN7_g4p3_advnoise_20250212.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/litec_filtered_CRN7_g4p3_advnoise_20250212.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-
#SBATCH --array=0-8
#SBATCH --requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

Niter=800000
gwb_Nfreqs=7
dataset=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/litec_filtered_ePSRs
outdir=/vast/palmer/scratch/mingarelli/bbl29/IPTA_DR2_analysis/chains/dr2litec_filtered/CRN${gwb_Nfreqs}_g4p3_advnoise
emp_dist_path=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/empdists/dr2full_CRN13_g4p3_advnoise.pkl
noisedict_path=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/noisedicts/dr2full_advnoise.json


srun -n $SLURM_NTASKS python3 ./scripts/run_gwb_adv_noise.py --fixed_gamma --gwb_Nfreqs $gwb_Nfreqs --dataset $dataset --Niter $Niter --outdir $outdir --emp_dist_path $emp_dist_path --noisedict_path $noisedict_path --human bjorn.larsen

