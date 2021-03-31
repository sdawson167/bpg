#!/bin/bash

#SBATCH --time 0-12:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3G
#SBATCH --output="01_testRanges04_sig.out"
#SBATCH --job-name="tau_04_01s"
#SBATCH -A def-shi
 
TABLE="sig_gamma_ranges"
N=$(cat $TABLE | wc -l)

tau=-0.4

for ((i=1; i<=N; i++))
do
	LINE=$(sed -n ${i}p $TABLE)
	read phase gamma0 gammaEnd <<< $LINE

	./BPG_TEST PW 3 $tau $gamma0 $gammaEnd 0.0 1.0 $phase  
done
