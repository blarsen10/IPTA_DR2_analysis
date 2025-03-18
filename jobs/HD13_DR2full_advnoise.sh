#!/bin/bash
#SBATCH --job-name=HD13_DR2full_advnoise_20241119
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/HD13_DR2full_advnoise_20241119.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/HD13_DR2full_advnoise_20241119.txt
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

gwb_Nfreqs=13
coredir=/home/bbl29/project/IPTA_DR2_analysis/dr2full/CRN${gwb_Nfreqs}_advnoise
dataset=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/full_ePSRs
noisedict=/home/bbl29/IPTA_DR2_analysis/noisedicts/dr2full_advnoise.json
Niter=10000

python3 ./scripts/run_reweighting.py --coredir $coredir --dataset $dataset --gwb_Nfreqs $gwb_Nfreqs --noisedict_path $noisedict -N $Niter

