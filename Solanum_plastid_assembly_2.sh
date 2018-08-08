#!/bin/sh
### This is step 2 in the plastid assembly. It takes the draft genome and looks for SNPs and InDels (from GATK) to look for mistakes during the assembly process. The output is used as an input into Solanum_plastid_assembly_3.sh.
# chmod 755 Solanum_plastid_assembly_2.sh
# nohup ./Solanum_plastid_assembly_2.sh &
###Note: Reference sequence only contains 1 inverted repeat

#####EDIT HERE#####
PI="amb365317"
start_of_IR="85942"
end_of_IR="111535"
start_of_IR2="129894"
final_length="155487"
#####EDIT DONE#####

##########OTHER PARAMETERS##########
###Length of LSC, IR, SSC;
LSC_length=$((start_of_IR - 1))
IR_length=$((end_of_IR - start_of_IR +1))
SSC_length=$(((start_of_IR2 - 1) - (end_of_IR +1) +1 ))
total_length=$((LSC_length + (2 * IR_length) + SSC_length))
###Check if lengths were done properly
if [[ $final_length -ne $total_length ]] ; then
    exit 1
fi
ref_plastid="${PI}.plastid.fa"
###Location of Fastq files
fastq1="/path/${PI}_1_clean.fq.gz"
fastq2="/path/${PI}_2_clean.fq.gz"
###RG header
SM=`echo $PI | cut -c 4-9`
ID=`echo $PI | cut -c 1-9`
Read_Group="@RG\tID:${ID}\tLB:barcode\tPL:illumina\tSM:${SM}\tPI:150"
##########END PARAMETERS##########


##############################START##############################
###Map all reads to the final plastid sequence (by accession) using bwa
if [ ! -s ${PI}.plastid_mpq.bam ]; then #check if already done
	bwa index ${ref_plastid}
	bwa mem -a -M -t 8 -R ${Read_Group} ${ref_plastid} ${fastq1} ${fastq2} >${PI}.plastid.sam
	###sort and index
	samtools faidx ${ref_plastid}
	samtools view -bt ${ref_plastid}.fai -o ${PI}.plastid.bam ${PI}.plastid.sam
	samtools sort ${PI}.plastid.bam ${PI}.plastid_s
	samtools view -b -q 30 -F 8 ${PI}.plastid_s.bam >${PI}.plastid_mpq.bam
	samtools index ${PI}.plastid_mpq.bam
	rm ${PI}.plastid.sam
	rm ${PI}.plastid.bam
	rm ${PI}.plastid_s.bam
else
	echo "plastid sequence has been mapped, sorted, and indexed"
fi

###Use GATK to confirm that the sequence is correct
if [ ! -s ${PI}.UGvariants.vcf ] ; then #check if done
	java -jar -Xmx4g /usr/local/bin/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=${ref_plastid} OUTPUT=${ref_plastid::(-2)}dict
	java -Xmx2g -jar /usr/local/bin/picard-tools/picard.jar MarkDuplicates I=${PI}.plastid_mpq.bam M=metrics.${PI}.txt TMP_DIR=/dauc2/tmp MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 O=${PI}.mpqu.bam
	samtools index ${PI}.mpqu.bam
	java -jar -Xmx8g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref_plastid} -I ${PI}.mpqu.bam -o ${PI}plastid.target_intervals.list
	java -jar -Xmx8g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref_plastid} -I ${PI}.mpqu.bam -targetIntervals ${PI}plastid.target_intervals.list -o ${PI}plastid.realigned.bam
	java -jar -Xmx4g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${ref_plastid} -I ${PI}plastid.realigned.bam -ploidy 2 -glm BOTH -baq CALCULATE_AS_NECESSARY -dt NONE -o ${PI}.UGvariants.vcf
else
	echo "GATK was already ran"
fi

###Use GATK to find coverage of plastid (LSC,1 IR, SSC)
if [ ! -s ${PI}_plastid_cov ]; then #check if done
	java -jar -Xmx4g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T DepthOfCoverage -R ${ref_plastid} -o ${PI}_plastid_cov -I ${PI}plastid.realigned.bam
else
	echo "Plastid coverage has been completed"
fi

###divide Plastid by region and find average coverage
###Then print out stats
if [ ! -s ${PI}_plastid_covLSC ]; then #check if done
	###separate file by region
	sed -n "2,$((LSC_length+1))p" ${PI}_plastid_cov >${PI}_plastid_covLSC
	sed -n "$((LSC_length+2)),$((LSC_length+2+IR_length))p" ${PI}_plastid_cov >${PI}_plastid_covIR
	sed -n "$((LSC_length+3+IR_length)),$((LSC_length+3+IR_length+SSC_length))p" ${PI}_plastid_cov >${PI}_plastid_covSSC
	
	###Find the average coverage of region
	echo "Average coverage of LSC is " >>${PI}plastid_final_stats
	awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}_plastid_covLSC >>${PI}plastid_final_stats
	echo "Average coverage of IR is " >>${PI}plastid_final_stats
	awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}_plastid_covIR >>${PI}plastid_final_stats
	echo "Average coverage of SSC is " >>${PI}plastid_final_stats
	awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}_plastid_covSSC >>${PI}plastid_final_stats

	###Print out regions that do not meet the Low (15% (LSC or SSC) or 30% mito (IR)) and High (1.5X (LSC or SSC) and 3X (IR)) regions
    Plastid_max=`awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}_plastid_covLSC`
    Pint=`printf "%.0f\n" $Plastid_max`
    MmaxSC=$((Pint * 15 / 100))
    PmaxSC=$((Pint * 3 / 2))
    MmaxIR=$((Pint * 3 / 10))
    PmaxIR=$((Pint * 3))
    awk '$2<'$MmaxSC ${PI}_plastid_covLSC >>${PI}_plastid_LowRegions
    awk '$2>'$PmaxSC ${PI}_plastid_covLSC >>${PI}_plastid_HighRegions
    awk '$2<'$MmaxSC ${PI}_plastid_covSSC >>${PI}_plastid_LowRegions
    awk '$2>'$PmaxSC ${PI}_plastid_covSSC >>${PI}_plastid_HighRegions
    awk '$2<'$MmaxIR ${PI}_plastid_covIR >>${PI}_plastid_LowRegions
    awk '$2>'$PmaxIR ${PI}_plastid_covIR >>${PI}_plastid_HighRegions
else
	echo "Plastid stats have been completed"
fi

###Find those variants that the alternate is greater than the reference
read -r -d '' perlscript_filter_plastid_vcf <<'EOF1'
use strict; use warnings;

my ($vcf) = @ARGV;

my (@col, @colon, @comma);

{
if ($vcf =~ /.gz$/) {open(FH, "zcat $vcf |") or die;}
else {open(FH,"<","$vcf") or die;}

while (my $lines = <FH>) {
  chomp $lines;
  next if ($lines =~ m/^#/); #Skip header lines with "#"
  my @col = split (/\t/, $lines); #Split the rest of the lines by tab
  my @colon = split (/:/, $col[9]); #Split the last tab by a :
  my @comma = split (/,/, $colon[1]); #Split the second variable by a ,
  if ($comma[1] > $comma[0]) {print "$lines\n"};
}
}
EOF1
perl -e "$perlscript_filter_plastid_vcf" ${PI}.UGvariants.vcf >${PI}.variants.vcf

###Create Coverage plot vs GC content (101 size bins)
bb.coverage --reference=${ref_plastid} --bam=${PI}plastid.realigned.bam --bed=${PI}.plastid.bed --gc --image=${PI}.plastid_CovGC

mv *CovGC.png ${PI}.CovGC.png

rm ${PI}_plastid_cov.sample_statistics ${PI}_plastid_cov.sample_summary ${PI}_plastid_cov.sample_cumulative_coverage_counts ${PI}_plastid_cov.sample_interval_summary ${PI}_plastid_cov.sample_interval_statistics ${PI}_plastid_cov.sample_cumulative_coverage_proportions ${PI}_plastid_covLSC ${PI}_plastid_covIR ${PI}_plastid_covSSC *.bed ${PI}.mpqu.bam* ${PI}.plastid_mpq.bam* metrics* *target_intervals.list *.fa.* *.dict *.data ${PI}.UGvariants.vcf*
