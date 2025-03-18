#!/bin/bash
#SBATCH --job-name=DR2lite_unfiltered_J1903+0327_factlike_20250303
#SBATCH --partition=scavenge
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/IPTA_DR2_analysis/logs/DR2lite_unfiltered_J1903+0327_factlike_20250303.txt
#SBATCH --error=/home/bbl29/IPTA_DR2_analysis/logs/error/DR2lite_unfiltered_J1903+0327_factlike_20250303.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu
#SBATCH --requeue

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

gwb_Nfreqs=13
n_live=1000
n_networks=4
projectpath=/home/bbl29/IPTA_DR2_analysis
outdir=/vast/palmer/home.grace/bbl29/project/IPTA_DR2_analysis/dr2lite_unfiltered/CRN${gwb_Nfreqs}_FL
noisepath=$projectpath/noisedicts/dr2lite_unfiltered_advnoise.json
psr=/vast/palmer/home.grace/bbl29/IPTA_DR2_analysis/data/lite_unfiltered_53_ePSRs/J1903+0327.pkl


python3 ./scripts/run_nautilus_factlike.py --gwb_Nfreqs $gwb_Nfreqs --Ncpus $SLURM_NTASKS --Ncpus $SLURM_NTASKS --Nlive $n_live --Nnet $n_networks --psr $psr --outdir $outdir --noisepath $noisepath

