#!/bin/bash
#
#SBATCH --partition=priority
#SBATCH --job-name=MIARNA
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=100000
#SBATCH --time=3-00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=benjamin.finander@childrens.harvard.edu

module load gcc/6.2.0
module load python/2.7.12

MASTER_DIR='/n/scratch3/users/b/bf78/MIARPJune15'
GENOME_DIR='/n/groups/neuroduo/lisa/reference_files/STAR_GRCm38_genome_dir'
annotation_gtf='/n/groups/neuroduo/Erin/reference_files.dir/Gencode_m19_GRCm38.p6/gencode.vM19.annotation.gtf'
TDF_GENOME='/n/groups/neuroduo/Erin/reference_files.dir/Gencode_m19_GRCm38.p6/GRCm38.p6.genome'

cd $MASTER_DIR

#replace with your sample codes from basespace
for sample in BK-MIA1 BK-MIA2 BK-MIA3 BK-MIA4 BK-MIA5 BK-MIA6 BK-MIA7 BK-MIA8 BK-PBS1 BK-PBS2 BK-PBS3 BK-PBS4 BK-PBS5 BK-PBS6 BK-PBS7 BK-PBS8 18PBS4M 18PBS1M 14PBS8F 14PBS6F 14PBS4M 14PBS2F 14PBS29M 14PBS28M 14PBS22M 14PBS10F 14MIA6M 14MIA4M 14MIA2F 14MIA25F 14MIA23F 14MIA11M 14MIA10M
do



# align trimmed reads to the  genome
module load star/2.5.4a
    STAR \
        --genomeDir $GENOME_DIR \
        --runThreadN 20 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --readFilesIn "$sample"_1.fastq "$sample"_2.fastq \
        --outFileNamePrefix "$sample" \
        --outSAMattributes NH HI NM MD \
        --outSAMtype BAM Unsorted

# Sort the reads
module load picard/2.8.0

    java -jar $PICARD/picard-2.8.0.jar SortSam \
        I="$sample"Aligned.out.bam \
        O="$sample"_sort.bam \
        SORT_ORDER=coordinate

# make tracks

    echo "track type=bedGraph name=\" $sample $1 "plus \" description=\" $sample "positive strand\" visibility=full autoScale=off windowingFunction=maximum viewLimits=0:30 color=150,0,75 color=0,0,150 priority=20" > "$sample".pos.bedGraph
    echo "track type=bedGraph name=\" $sample $1 "minus \" description=\" $sample "minus strand\" visibility=full autoScale=off negateValues=on windowingFunction=maximum viewLimits=-30:0 color=150,0,75 color=0,0,150 priority=20" > "$sample".min.bedGraph

    STAR \
        --runMode inputAlignmentsFromBAM \
        --inputBAMfile "$sample"_sort.bam \
        --outWigType bedGraph \
        --outWigNorm None \
        --outWigStrand Stranded

    cat Signal.Unique.str1.out.bg >> "$sample".min.bedGraph
    cat Signal.Unique.str2.out.bg >> "$sample".pos.bedGraph

#And make tdf files (a binary file that loads easily into UCSC and IGV, similar to .bw but smaller)
	module load igvtools/2.3.98

        java -jar $IGVTOOLS/igvtools.jar toTDF \
		-f mean,max \
		"$sample".pos.bedGraph \
		"$sample".pos.tdf \
		$TDF_GENOME

	java -jar $IGVTOOLS/igvtools.jar toTDF \
		-f mean,max \
		"$sample".min.bedGraph \
		"$sample".min.tdf \
		$TDF_GENOME

# Make the master counts:
#When using the Clontech RNA-seq library prep sequenced on NextSeq, counts should be on reverse, not forward strand
module load htseq/0.9.1

#I like to quantify both the gene-level and exon-level counts for s4U metabolic labeling (s4U-RNA should be enriched for intronic reads). You probably just need the exon-level counts

    htseq-count -f bam \
        -t exon \
        -i gene_name \
        -s reverse \
	--minaqual=5 \
        "$sample"_sort.bam \
        $annotation_gtf > mature."$sample".htseq.txt

done


