#!/bin/bash
#SBATCH --job-name=DR2full_J1022+1001_factlike_20241101
#SBATCH --partition=scavenge
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/DR2full_J1022+1001_factlike_20241101.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/DR2full_J1022+1001_factlike_20241101.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

n_live=1000
n_networks=4
projectpath=/home/bbl29/IPTA_DR2_analysis
outdir=/vast/palmer/home.grace/bbl29/project/IPTA_DR2_analysis/dr2full/factlike
noisepath=$projectpath/noisedicts/dr2full_advnoise.json
psr=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/full_ePSRs/J1022+1001.hdf5


python3 ./scripts/run_nautilus_factlike.py --Ncpus $SLURM_NTASKS --Ncpus $SLURM_NTASKS --Nlive $n_live --Nnet $n_networks --psr $psr --outdir $outdir --noisepath $noisepath

