#!/bin/bash
#SBATCH -J MTS_size 		# Job name
#SBATCH -o MTS_size.o%j 	# Name of stdout output file (%j expands to jobId)
#SBATCH -e MTS_size.o%j 	# Name of stderr output file(%j expands to jobId)
#SBATCH -n 128 		   		# Total number of mpi tasks requested
#SBATCH -p development	 	# Submit to the 'normal' or 'development' queue
#SBATCH -t 1:00:00 			# Run time (hh:mm:ss)
#SBATCH --mail-user=jepusto@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < "MTS sample size sims - TACC".R 
