#!/usr/bin/env bash
# Run ft_load

#SBATCH --output=ft_load-%j.out
#SBATCH -t 30:00
#SBATCH --nodes=4
#SBATCH --mem-per-cpu 25000

# Load modules
module load MATLAB/2018b

# Go to directory
cd /gpfs/milgram/project/turk-browne/projects/stimulation_behavior/aal66stim/stimulation_network/ERPanalysis

# Run script
matlab -nodesktop -nojvm -r "addpath(pwd); addpath('/gpfs/milgram/project/turk-browne/projects/stimulation_behavior/intermediate_data'); dataERPfix; exit"
