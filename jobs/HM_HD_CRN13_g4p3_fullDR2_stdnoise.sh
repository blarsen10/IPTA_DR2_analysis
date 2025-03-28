#!/bin/bash
#SBATCH --job-name=HM_GWB13_g4p3_stdnoise_20240719
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/HM_GWB13_g4p3_stdnoise_20240719.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/HM_GWB13_g4p3_stdnoise_20240719.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=3-
#SBATCH --array=0-7
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

Niter=800000
dataset=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/full_ePSRs
outdir=/vast/palmer/scratch/mingarelli/bbl29/IPTA_DR2_analysis/chains/dr2full/HM_GWB13_g4p3_stdnoise
emp_dist_path=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/empdists/dr2full_crn_stdnoise.pkl
noisedict_path=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/noisedicts/dr2full_stdnoise.json


srun -n $SLURM_NTASKS python3 ./scripts/run_gwb_hm.py --fixed_gamma --dataset $dataset --Niter $Niter --outdir $outdir --emp_dist_path $emp_dist_path --noisedict_path $noisedict_path --human bjorn.larsen

