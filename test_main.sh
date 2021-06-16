#!/bin/bash

#SBATCH --time 0-12:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3G
#SBATCH --output="test_tau_-0.1_sig.out"
#SBATCH --job-name="test_tau_-0.1_sig"
#SBATCH -A def-shi
 
tau=-0.1

./BPG_TEST $tau 7  
