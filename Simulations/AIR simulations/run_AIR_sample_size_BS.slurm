#!/bin/bash
#SBATCH -J AIR_BS 		# Job name
#SBATCH -o AIR_BS.o%j 	# Name of stdout output file (%j expands to jobId)
#SBATCH -e AIR_BS.o%j 	# Name of stderr output file(%j expands to jobId)
#SBATCH -n 128		   		# Total number of mpi tasks requested
#SBATCH -p development	 	# Submit to the 'normal' or 'development' queue
#SBATCH -t 1:00:00 			# Run time (hh:mm:ss)
#SBATCH --mail-user=dms031000@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < "AIR bootstrap performance sims - TACC".R 
