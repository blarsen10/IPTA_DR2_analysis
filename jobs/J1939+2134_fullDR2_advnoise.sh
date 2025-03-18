#!/bin/bash
#SBATCH --job-name=DR2_J1939+2134_advnoise_20240925
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/J1939+2134_DR2full_advnoise_20240925.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/J1939+2134_DR2full_advnoise_20240925.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-
#SBATCH --array=0-7
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu
#SBATCH --requeue

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

psrname=J1939+2134
Niter=350000
outdir=/vast/palmer/scratch/mingarelli/bbl29/IPTA_DR2_analysis/chains/dr2full/advnoise/J1939+2134
psr=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/full_ePSRs/J1939+2134.hdf5
emp_dist_path=


srun -n $SLURM_NTASKS python3 ./scripts/run_custom_noise.py --psrname $psrname --Niter $Niter --outdir $outdir --psr $psr --emp_dist_path $emp_dist_path --human bjorn.larsen

