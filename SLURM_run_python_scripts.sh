#!/bin/sh
#SBATCH --account=def-wilsyl
#SBATCH --ntasks=1               # number of MPI processes
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G      # memory; default unit is megabytes
#SBATCH --time=0-06:00           # time (DD-HH:MM)

python3 extract_spatial_tuning.py
python3 select_neurons.py
#python3 seqNMF_extract_seq_results.py
#python3 seqNMF_extract_seqReplay_results.py
#python3 seqNMF_optimize_L.py
