#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=1-00:15:00     # 1 day and 15 minutes
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu
MRO_DISK_SPACE_CHECK=disable

export PATH=/rhome/rrugg002/bigdata/cellranger-6.1.2:$PATH

cellranger multi --id=WT_Female --csv=Relma_WT_Female_HFD.csv
