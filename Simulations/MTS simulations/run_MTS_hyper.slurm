#!/bin/bash
#SBATCH -J MTS_hyp 		# Job name
#SBATCH -o MTS_hyp.o%j 	# Name of stdout output file (%j expands to jobId)
#SBATCH -e MTS_hyp.o%j 	# Name of stderr output file(%j expands to jobId)
#SBATCH -n 64    		# Total number of mpi tasks requested
#SBATCH -p normal	 	# Submit to the 'normal' or 'development' queue
#SBATCH -t 3:00:00 		# Run time (hh:mm:ss)
#SBATCH --mail-user=jepusto@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < "MTS hyperparameter sims - TACC".R 
