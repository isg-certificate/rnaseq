#!/bin/bash
#SBATCH --job-name=get_genome
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Download genome and annotation from ENSEMBL
#################################################################

# load software
module load samtools/1.16.1

# output directory
GENOMEDIR=../../genome
mkdir -p $GENOMEDIR

# we're using Fundulus heteroclitus from ensembl v112
    # we'll download the genome, GTF annotation and transcript fasta
    # https://useast.ensembl.org/Fundulus_heteroclitus/Info/Index

# note that in these URLs we are downloading v112 specifically. 

# download the genome
wget -P ${GENOMEDIR} http://ftp.ensembl.org/pub/release-112/fasta/fundulus_heteroclitus/dna/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.dna_sm.toplevel.fa.gz

# download the GTF annotation
wget -P ${GENOMEDIR} http://ftp.ensembl.org/pub/release-112/gtf/fundulus_heteroclitus/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gtf.gz

# download the transcript fasta
wget -P ${GENOMEDIR} http://ftp.ensembl.org/pub/release-112/fasta/fundulus_heteroclitus/cdna/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cdna.all.fa.gz

# decompress files
gunzip ${GENOMEDIR}/*gz

# generate simple samtools fai indexes 
samtools faidx ${GENOMEDIR}/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.dna_sm.toplevel.fa
samtools faidx ${GENOMEDIR}/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cdna.all.fa

