#!/bin/sh
#SBATCH --job-name=Devon_plot
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1

#SBATCH --ntasks-per-node=12
#SBATCH --time=00:20:00


. ~/miniforge3/etc/profile.d/conda.sh
conda activate plotting

ASSEMBLY=Devon_UK
RMOUTDIR=/lustre/scratch/daray/dRet_curation/bed_files

WORKDIR=/lustre/scratch/daray/dRet_curation/plotting_redo/${ASSEMBLY}_RM
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 

cd $WORKDIR

python $CURATIONPATH/te_plotting2.py \
	-r ${RMOUTDIR}/${ASSEMBLY}.fa.out.gz \
	-s ${RMOUTDIR}/${ASSEMBLY}.summary.gz \
	-op ${ASSEMBLY} \
	-c LINE,SINE,LTR,DNA,RC,Unknown \
	-proc 12
	
