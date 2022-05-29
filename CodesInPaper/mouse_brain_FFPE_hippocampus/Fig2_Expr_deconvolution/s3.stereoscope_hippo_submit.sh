#!/bin/bash
#SBATCH --job-name stereoscope
#SBATCH -p super
#SBATCH -N 1
#SBATCH --mem 252928
#SBATCH -t 0-100:0:00
#SBATCH -o /archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/res/log/stereoscope.%j.out
#SBATCH -e /archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/res/log/stereoscope.%j.err

module load python/3.8.x-anaconda
source activate myVE.python3.8
module load gcc/6.3.0

cd /archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/res

stereoscope run \
--sc_cnt ../data/MouseBrainHippo/sc_cnt.tsv \
--sc_labels ../data/MouseBrainHippo/sc_mta.tsv \
--st_cnt ../data/MouseBrainHippo/st_cnt.tsv \
-o MouseBrainHippo


