#!/bin/bash
#SBATCH -J MTS_BS 		# Job name
#SBATCH -o MTS_BS.o%j 	# Name of stdout output file (%j expands to jobId)
#SBATCH -e MTS_BS.o%j 	# Name of stderr output file(%j expands to jobId)
#SBATCH -n 64 		   		# Total number of mpi tasks requested
#SBATCH -p development	 	# Submit to the 'normal' or 'development' queue
#SBATCH -t 0:30:00 			# Run time (hh:mm:ss)
#SBATCH --mail-user=dms031000@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < "MTS bootstrap performance sims - TACC".R 
