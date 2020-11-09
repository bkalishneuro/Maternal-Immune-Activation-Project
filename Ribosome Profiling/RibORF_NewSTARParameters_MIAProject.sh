#!/bin/bash
#
#SBATCH --partition=priority
#SBATCH --job-name=Ribo_Hs
#SBATCH -c 8 #Using 8 cores
#SBATCH -N 1 #Using 1 node, and all cores on same node
#SBATCH --mem=80000
#SBATCH --time=1-00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=benjamin.finander@childrens.harvard.edu

#usage: sbatch ED_Hs_Riboseq.sh samplename

MASTER_DIR='/n/scratch3/users/b/bf78/LastMIARP'
fastq_dir1='/n/groups/neuroduo/ARCHIVED_neuroduo_Tier3/Erin/Ribosome_profiling.dir/Human_adult.dir/Ribo_seq.dir/200309_Broad_RP.dir/fastq'
fastq_dir2='/n/groups/neuroduo/ARCHIVED_neuroduo_Tier3/Erin/Ribosome_profiling.dir/Human_adult.dir/Ribo_seq.dir/200505_Broad_RP.dir/fastq'
ref_dir="/n/groups/neuroduo/Erin/reference_files.dir/Gencode_m19_GRCm38.p6"
STAR_gencode='/n/groups/neuroduo/brian/ben/ReferenceGenomes/STAR_GRCm38_genome_SpecialRiboseq'
STAR_lncRNAKB='/n/groups/neuroduo/Erin/reference_files.dir/grch38/STAR_lncRNAKB.dir'
RibORF="/n/groups/neuroduo/Erin/SOFT/RibORF-master/RibORF.1.0"
annotation_gtf='/n/groups/neuroduo/Erin/reference_files.dir/Gencode_m19_GRCm38.p6/gencode.vM19.annotation.gtf'
sample=$1
ref=$2

echo "$sample" "$ref"

mkdir $MASTER_DIR/RibORF_"$ref".dir/"$sample".dir
cp "$sample"_Ribo.fastq $MASTER_DIR/RibORF_"$ref".dir/"$sample".dir
cd $MASTER_DIR/RibORF_"$ref".dir/"$sample".dir

if [ $ref == "gencode" ]
then
        GENOME_DIR=$STAR_gencode
        RibORF_ref=$ref_dir/gencode.vM19.annotation.genePred
        RibORF_ORFs=$ref_dir/candidateORF.genepred.txt
else
        GENOME_DIR=$STAR_lncRNAKB
        RibORF_ref=$ref_dir/lncRNAKB_hg38_v7.genePred.txt
        RibORF_ORFs=$ref_dir/lncRNAKB_hg38.candidateORF.genepred.txt
fi

#Handles merged files earlier, which is probably better than merging after alignment.



#Trim with cutadapt - this is necessary for the Clontech smRNA kit because
#a poly-A tail is added before the RT and needs to be removed
module load gcc/6.2.0 python/2.7.12 cutadapt/1.14 R/3.5.1 bowtie2/2.3.4.3 star/2.5.4a

cutadapt \
    -m 15 \
    -a AAAAAAAAAA \
    -o trim."$sample".fastq \
    "$sample"_Ribo.fastq \

# Filter all rRNAs
module load bowtie2/2.3.4.3

bowtie2 \
    -x /n/groups/neuroduo/Erin/reference_files.dir/Gencode_m19_GRCm38.p6/mm10.ribosome \
    -U trim."$sample".fastq \
    --un norrna."$sample".fastq \
    -S ribosome."$sample".fastq

# Align non-rRNA reads to genome
#Extra parameters included for short read mapping (reference here: https://sites.google.com/site/dvanichkina/ngs#TOC-Working-with-very-short-reads-25-nt-as-for-example-from-ribosomal-profiling-)
STAR \
        --genomeDir $GENOME_DIR \
        --runThreadN 20 \
        --alignEndsType EndToEnd \
        --clip5pNbases 3 \
        --seedSearchStartLmax 15 \
        --outSJfilterOverhangMin 30 8 8 8 \
        --outFilterScoreMin 0 \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterMatchNmin 0 \
        --outFilterMatchNminOverLread 0.66 \
        --readFilesIn norrna."$sample".fastq \
        --outFileNamePrefix "$sample"_"$ref" \
        --outSAMattributes All \
        --outSAMtype BAM Unsorted

# Sort the reads
module load picard/2.8.0

java -jar $PICARD/picard-2.8.0.jar SortSam \
    I="$sample"_"$ref"Aligned.out.bam \
    O=Ribo."$sample"_"$ref"_sort.bam \
    SORT_ORDER=coordinate

rm "$sample"_"$ref"Aligned.out.bam

module load samtools/1.9
samtools view -h -o "$sample".sam Ribo."$sample"_"$ref"_sort.bam

#Make tracks
module load star/2.5.4a

STAR \
    --runMode inputAlignmentsFromBAM \
    --inputBAMfile Ribo."$sample"_"$ref"_sort.bam \
    --outWigType bedGraph \
    --outWigNorm None \
    --outWigStrand Stranded

cat Signal.Unique.str1.out.bg > Ribo."$sample"_"$ref".pos.bedGraph
cat Signal.Unique.str2.out.bg > Ribo."$sample"_"$ref".min.bedGraph

rm *.bg 

#Run RibORF on the pipeline
perl $RibORF/readDist.pl \
        -f "$sample".sam \
        -g $RibORF_ref \
        -o . \
        -d 19,20,21,22,28,29,30, \
       -l 40 \
       -r 70

#Make plots from read distribution output
R CMD BATCH readDist.plot.19.R
R CMD BATCH readDist.plot.20.R
R CMD BATCH readDist.plot.21.R
R CMD BATCH readDist.plot.22.R
R CMD BATCH readDist.plot.28.R
R CMD BATCH readDist.plot.29.R
R CMD BATCH readDist.plot.30.R

#####-----IMPORTANT: Be sure to check the periodicity before continuing. 
# Determine which read lengths have good periodicity (usually 28-30nt) and
# proper offset correction (usually 12). Create a file in $sample.dir called
#"offset.correction.parameters.txt which usually takes the format:
#21 12
#28 12
#29 12
#30 12

#Correct read locations based on offset distances between 5' ends and ribosomal A-sites
perl $RibORF/offsetCorrect.pl \
       -r $MASTER_DIR/RibORF_"$ref".dir/"$sample".dir/"$sample".sam \
        -p $MASTER_DIR/offset.correction.parameters.txt \
        -o corrected."$sample".sam

#And calculate read distribution from offset-corrected file
perl $RibORF/readDist.pl \
       -f corrected."$sample".sam \
        -g $RibORF_ref \
        -o $MASTER_DIR/RibORF_"$ref".dir/"$sample".dir \
        -d 1

#Run RibORF to identify translated ORFs
perl $RibORF/ribORF.pl \
        -f corrected."$sample".sam \
        -c $RibORF_ORFs \
        -o $MASTER_DIR/RibORF_"$ref".dir/"$sample".dir \
        -l 24

R CMD BATCH ribORF.learning.R

#The file you want with significantly translated ORFs and their counts is called
#repre.valid.pred.pvalue.parameters.txt
mv repre.valid.pred.pvalue.parameters.txt "$sample"."$ref".repre.valid.pred.pvalue.parameters.txt

module load htseq/0.9.1

#I like to quantify both the gene-level and exon-level counts for s4U metabolic labeling (s4U-RNA should be enriched for intronic reads). You probably just need the exon-level counts

    htseq-count -f bam \
        -t exon \
        -i gene_name \
    --minaqual=5 \
        Ribo."$sample"_"$ref"_sort.bam \
        $annotation_gtf > Ribo."$sample".htseq.txt