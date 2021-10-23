#!/bin/bash
#SBATCH --account=def-aschmidt  # replace this with your supervisors account
#SBATCH --nodes=1                # number of node MUST be 1
#SBATCH --cpus-per-task=4        # number of processes
#SBATCH --mem-per-cpu=8192M      # memory; default unit is megabytes
#SBATCH --time=3-00:00           # time (DD-HH:MM)
#SBATCH --mail-user=widemberg@dme.ufrj.br # Send email updates to you or someone else
#SBATCH --mail-type=ALL          # send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0
module load r/4.0.2
ml gcc r-bundle-bioconductor

# Export the nodes names. 
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export R_LIBS=~/local/R_libs/
R -f TBsetting_case2_Scen1_gera_X.R
