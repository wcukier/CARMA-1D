#!/bin/bash
#SBATCH --job-name=carma_test
#SBATCH --error=carma_test.err
#SBATCH --partition=caslake
#SBATCH --account=pi-diana8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wcukier.notify@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01-12:00:00

pwd; hostname; date

echo "Running program on $SLURM_JOB_NUM_NODES nodes, with each node getting $SLURM_NTASKS_PER_NODE tasks, for $SLURM_NTASKS total tasks."

module load intel

srun --exclusive ./run-carma.csh diamondback_test.exe > diamondback_test.out &

wait

date