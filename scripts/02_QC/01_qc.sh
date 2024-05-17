#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastqc/0.12.1
module load parallel/20180122

# set input/output directory variables
INDIR=../../data/fastq/
REPORTDIR=../../results/02_qc/fastqc_reports
mkdir -p $REPORTDIR

ACCLIST=../../metadata/accessionlist.txt
# run fastp in parallel, 10 samples at a time
cat $ACCLIST | parallel -j 10 \
    "fastqc --outdir $REPORTDIR $INDIR/{}_{1..2}.fastq.gz"