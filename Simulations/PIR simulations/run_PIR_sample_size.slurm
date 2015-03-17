#!/bin/bash
#SBATCH -J PIR_size 		# Job name
#SBATCH -o PIR_size.o%j 	# Name of stdout output file (%j expands to jobId)
#SBATCH -e PIR_size.o%j 	# Name of stderr output file(%j expands to jobId)
#SBATCH -n 128		   		# Total number of mpi tasks requested
#SBATCH -p normal		 	# Submit to the 'normal' or 'development' queue
#SBATCH -t 3:00:00 			# Run time (hh:mm:ss)
#SBATCH --mail-user=jepusto@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < "PIR sample size sims - TACC.R"
