#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --output=snakemake.log
#SBATCH --error=snakemake.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=node2

# (1) Conda 초기화
eval "$(conda shell.bash hook)"
conda activate snakemake


cd /data/Storage_DAS03/home3/menteeD/CWS/RNA/KBPMA_RNAseq-main


# Snakemake 실행
snakemake --cores 30 --use-conda --rerun-incomplete --printshellcmds

