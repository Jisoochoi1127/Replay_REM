#!/bin/bash
#SBATCH --job-name=Replay_batch
#SBATCH --account=def-wilsyl
#SBATCH --time=22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=124G
#SBATCH --mail-user=jisuchoi1127@gmail.com
#SBATCH --mail-type=ALL
module load StdEnv/2020
module load matlab/2020a
matlab -nodisplay -batch "MASTER_replay_REM_jump30('Selection')"
