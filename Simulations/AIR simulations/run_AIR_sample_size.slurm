#!/bin/bash
#SBATCH -J AIR_size 		# Job name
#SBATCH -o AIR_size.o%j 	# Name of stdout output file (%j expands to jobId)
#SBATCH -e AIR_size.o%j 	# Name of stderr output file(%j expands to jobId)
#SBATCH -n 64 		   		# Total number of mpi tasks requested
#SBATCH -p development	 	# Submit to the 'normal' or 'development' queue
#SBATCH -t 1:00:00 			# Run time (hh:mm:ss)
#SBATCH --mail-user=jepusto@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < "AIR sample size sims - TACC".R 
