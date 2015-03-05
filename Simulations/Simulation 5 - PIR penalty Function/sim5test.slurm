#!/bin/bash
#SBATCH -J sim5test
#SBATCH -o sim5test.o%j
#SBATCH -e sim5test.o%j
#SBATCH -n 32
#SBATCH -p development
#SBATCH -t 01:00:00
#SBATCH -A Parallel-demonstrati
#SBATCH --mail-user=dms031000@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats
ibrun RMPISNOW < sim5test.R
