#!/bin/bash
#$ -N script_fastqc_raw_menta
#$ -cwd
#$ -V
#$ -l h_vmem=3G
#$ -S /bin/bash
#$ -pe PE 1

cd /home/llopez/GATKkorea/reference
# Genome indexing for BWA
bwa index /home/llopez/GATKkorea/reference/Pvulgaris_442_v2.0.fa
# Creation of 2 new index files for reference (in formats used by Samtools and GATK)
samtools faidx /home/llopez/GATKkorea/reference/Pvulgaris_442_v2.0.fa
java -jar -Xmx10G  /home/llopez/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar CreateSequenceDictionary -R /home/llopez/GATKkorea/reference/Pvulgaris_442_v2.0.fa -O /home/llopez/GATKkorea/reference/Pvulgaris_442_v2.0.dict


## Building g.VCF files of each sample per GATK4.0 protocol
# alignment by bwa

for i in /home/llopez/trimming/KOREA/Trimmomatic/Trimmomatic_on_*_1.fastq
do
	# Locating the destination of the process
	cd /home/llopez/GATKkorea
	SUBSTRING=$(echo $i| cut -d'_' -f 3)
	echo `date`
	echo "finished upload"
	# alignment by bwa
	## example shows: $SUBSTRING.fastq.gz
	bwa mem -R '@RG\tID:$SUBSTRING\tSM:$SUBSTRING\tPL:Illumina' /home/llopez/GATKkorea/reference/Pvulgaris_442_v2.0.fa \
	$i > $SUBSTRING.sam
	#clean sam files
	java -jar -Xmx10G  /home/llopez/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar CleanSam -INPUT $SUBSTRING.sam -VALIDATION_STRINGENCY SILENT -OUTPUT $SUBSTRING.clean.sam
	rm $SUBSTRING.sam
	#This tool ensures that all mate information is in sync between each read and its mate.
	java -jar -Xmx10G  /home/llopez/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar FixMateInformation -INPUT $SUBSTRING.clean.sam -VALIDATION_STRINGENCY SILENT -OUTPUT $SUBSTRING.clean.fix.sam
	rm $SUBSTRING.clean.sam
	#Convert from sam to bam
	samtools view -S $SUBSTRING.clean.fix.sam -b -o $SUBSTRING.clean.fix.bam
	rm $SUBSTRING.clean.fix.sam
	#order alignment
	#via samtools
	#samtools sort -o G1_1.sort.bam -O BAM  G1_1.bam
	#via gatk
	java -jar -Xmx10G  /home/llopez/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SortSam -INPUT $SUBSTRING.clean.fix.bam -OUTPUT $SUBSTRING.sort.bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY SILENT
	rm $SUBSTRING.clean.fix.bam
	#statistics
	#chmod 0777 Bamtools/1.0
	#bamtools stats -in $SUBSTRING.sort.bam > $SUBSTRING.sort.bam.STATS_CON_DUP
	#Remove duplicates
	java -jar -Xmx10G  /home/llopez/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar MarkDuplicates -I $SUBSTRING.sort.bam -O $SUBSTRING.dup.bam -REMOVE_DUPLICATES TRUE  -METRICS_FILE $SUBSTRING.dup.metrics.txt
	rm $SUBSTRING.sort.bam
	#Statistics after cleaning
	bamtools stats -in $SUBSTRING.dup.bam > $SUBSTRING.dup.bam.STATS_SIN_DUP
	# index for gatk
	samtools index $SUBSTRING.dup.bam
	# second file is dictionary reference
	#via samtools
	#rm GCF_000499845.1_PhaVulg1_0_genomic.dict
	#samtools dict -o /home/llopez/GATK/reference/GCA_001642375.1_Mlong1.0_genomic.dict /home/llopez/GATK/reference/GCA_001642375.1_Mlong1.0_genomic.fa
	#via gatk
	#rm GCF_000499845.1_PhaVulg1_0_genomic.dict
	# Get variants (gvcf) for a sample
	java -jar -Xmx10G  /home/llopez/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar HaplotypeCaller -I $SUBSTRING.dup.bam -R /home/llopez/GATKkorea/reference/Pvulgaris_442_v2.0.fa -O $SUBSTRING.g.vcf -ERC GVCF -ploidy 2
	rm $SUBSTRING.dup.bam
	echo `date`
	echo "finished obtaining g.VCF of each sample"
	#rm $SUBSTRING.g.vcf
done

# After the loop in each sample continue to collapse the gvcf to vcf
java -jar -Xmx10G  /home/llopez/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenotypeGVCFs -R /home/llopez/GATKkorea/reference/Pvulgaris_442_v2.0.fa -V Sample1.g.vcf -V Sample2.g.vcf -V Sample86.g.vcf -O All_samples.vcf
echo `date`
echo "completed obtaining molecular variants"
