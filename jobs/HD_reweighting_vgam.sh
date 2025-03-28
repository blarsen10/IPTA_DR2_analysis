#!/bin/bash
#SBATCH --job-name=reweight_GWB13_stdnoise_20240912
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/reweight_GWB13_stdnoise_20240912.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/reweight_GWB13_stdnoise_20240912.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-
#SBATCH --array=0-15
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

coredir=/home/bbl29/project/IPTA_DR2_analysis/dr2full/CRN13_stdnoise
dataset=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/full_ePSRs
noisedict=/home/bbl29/IPTA_DR2_analysis/noisedicts/dr2full_stdnoise.json
Niter=10000

python3 ./scripts/run_reweighting.py --coredir $coredir --dataset $dataset --noisedict_path $noisedict -N $Niter --fixed_gamma