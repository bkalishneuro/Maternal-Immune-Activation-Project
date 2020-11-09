#!/bin/bash
#
#SBATCH --partition=short
#SBATCH --job-name=REDO4
#SBATCH -c 4
#SBATCH --mem-per-cpu=20000
#SBATCH --time=0-05:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=benjamin.finander@childrens.harvard.edu

module load gcc/6.2.0
module load bcl2fastq/2.17.1.14

#For each sequencing run, change to the main directory, and run bcl2fastq to convert .bcl files to .fastq files

SequencingRunFolder=x
cd ${SequencingRunFolder} #Make sure you cd to the directory with the RunInfo.xml file, it will find the BCL files buried in there, don't worry
bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0