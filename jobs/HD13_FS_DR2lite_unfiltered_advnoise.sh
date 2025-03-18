#!/bin/bash
#SBATCH --job-name=HD13_FS_DR2lite_uf_43_advnoise_20250216
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/HD13_FS_DR2lite_uf_43_advnoise_20250216.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/HD13_FS_DR2lite_uf_43_advnoise_20250216.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-
#SBATCH --array=0-7
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

gwb_Nfreqs=13
coredir=/home/bbl29/project/IPTA_DR2_analysis/dr2lite_unfiltered_43/CRN${gwb_Nfreqs}_FS_advnoise
dataset=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/lite_unfiltered_43_ePSRs
noisedict=/home/bbl29/IPTA_DR2_analysis/noisedicts/dr2lite_unfiltered_advnoise.json
Niter=20000

python3 ./scripts/run_reweighting.py --coredir $coredir --dataset $dataset --gwb_Nfreqs $gwb_Nfreqs --noisedict_path $noisedict -N $Niter --freespec

