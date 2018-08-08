#!/bin/sh
### This is step 2 in the plastid assembly. This takes the output from Solanum_plastid_assembly_2.sh and automates the correction to the plastid sequences. Then, it checks for errors with GATK (again) and softclips. Any unusual softclips are manually corrected. If there are any remaining SNPs and small InDels (from GATK) or any softclips manually corrected, this script needs to be re-ran until everything comes back good.
# chmod 755 Solanum_plastid_assembly_3.sh
# nohup ./Solanum_plastid_assembly_3.sh &
###Note: Reference sequence only contains 1 inverted repeat

###List of accessions here
accessions="aba458403 aba458404 ach558032 acr365313 acr365314 acr498204 alb498206 amb365317 amb365362 amb498209 amb498210 amb498212 amb498213 and320345 and561648 and561658 avi498091 avi498092 avi498093 ber218225 ber498105 ber527886 ber545850 bla498214 bra320265 bra558460 bre310931 bre473378 bre498111 bre498218 bre545968 bre545970 bre545971 bre545981 buk266385 buk365353 buk414155 buk473492 buk473493 buk473494 buk568933 buk568954 bul545751 bul604074 caj230522 can210035 can246533 can265864 can265865 can283084 can442696 can473355 can498226 can498227 can545972 can568969 car283062 car283063 car347759 cha275138 cha320294 cha472816 cha472830 cha500020 cho365328 cho365339 etu498311 gon195186 gon195188 gon195214 gon458393 gou472911 gou472991 gou472995 gou473019 gou473077 gou473106 gou500022 gou537026 gou545865 gou545975 gou545978 gou558067 haw320364 hon473365 hon498067 hon498071 hon545879 hua320327 hua498255 hyp473477 inc473060 inc473067 inc473069 inc473070 inc500048 jam641944 jam664024 kur472924 kur472936 kur472948 kur472952 kur558185 kur558208 lax283088 lax498252 lax607887 lep458378 lep473342 lep473446 lep473451 lep545985 lep545987 lim473468 mar210040 mar310944 med210045 med230507 med320260 med458402 med473496 meg210034 meg473158 meg500029 meg546000 mul210044 mul210052 mul210055 mul275272 mul365336 mul365337 mul365338 mul473349 mul473352 mul498266 pal245763 pam275274 pam275275 pam442697 pam458381 pau473489 phu195191 phu195198 phu225665 phu225693 phu225703 phu243467 phu243468 phu243469 phu258855 pin253214 pin537023 pol161728 pol347770 sog230510 sog365360 spa246536 spa473375 spa473385 spa498134 spa498284 spa498285 spe320299 spe458335 spe458337 spe472966 spe472988 spe472990 ste195204 ste205527 ste230512 ste230513 ste234011 ste255527 ste283141 ste365344 tar217457 tar414152 tar458366 tar473217 tar473218 ver195170 ver275256 ver275260 ver320332 ver320333 ver458370 ver473303 ver473309 ver498010 ver498061 ver500070 ver545745 ver545747 ver545884 ver558150 ver558463 ver558488 vio473396 vio473398 vio498296"

for PI in $accessions ; do ###Loop through the list of accessions
mkdir ${PI}

##########OTHER PARAMETERS##########
ori_ref="/path/S.tuberosum_chloroplast.fa"
ref_plastid="${PI}.plastid2.fa"
###Location of Fastq files
fastq1="/path/${PI}_1_clean.fq.gz"
fastq2="/path/${PI}_2_clean.fq.gz"
###RG header
SM=`echo $PI | cut -c 4-9`
ID=`echo $PI | cut -c 1-9`
Read_Group="@RG\tID:${ID}\tLB:barcode\tPL:illumina\tSM:${SM}\tPI:150"
##########END PARAMETERS##########


##############################START##############################
###Fix the fasta file using the variant.vcf file
read -r -d '' perlscript_CorrectPlastid_from_vcf <<'EOF1'
use strict; use warnings;
my ($fasta, $vcf) = @ARGV;

open FILE, $fasta or die;
my $hdr = <FILE>;
chomp $hdr;
my $sequence = <FILE>;
chomp $sequence;
close FILE;

open FILE2, $vcf or die;
my @lines = reverse <FILE2>;
foreach my $line (@lines) {
    my @col = split (/\t/, $line);
    my $reflength = length($col[3]);
    my $position = $col[1] -1;
    if ($col[4] !~ m/,/) {substr($sequence,$position,$reflength)="$col[4]";}
    else {
        my @altbase = split (/,/,$col[4]);
        my @AD = split (/:/,$col[9]);
        my @altcoverage = split (/,/,$AD[1]);
        if ($altcoverage[2] > $altcoverage[1]) {substr($sequence,$position,$reflength)="$altbase[1]";}
        else {substr($sequence,$position,$reflength)="$altbase[0]";}
    }
}
print "$hdr\n$sequence\n";
EOF1

perl -e "$perlscript_CorrectPlastid_from_vcf" /dauc4/Binquan/plastid/plastid_first_check/${PI}/${PI}.plastid.fa /dauc4/Binquan/plastid/plastid_first_check/${PI}/${PI}.variants.vcf >${PI}/${PI}.plastid2.fa

###Find new coordinates of LSC, IR, and SSC
read -r -d '' perlscript_IR_sizes <<'EOF2'
use strict; use warnings;

my ($sizes, $changes) = @ARGV;

open FILE, $sizes or die;
my ($start, $end, $length);
while (<FILE>) {
    if (m/start_of_IR="(\d+)"/) {$start = $1;}
    elsif (m/end_of_IR="(\d+)"/) {$end = $1;}
    elsif (m/start_of_IR2="(\d+)"/) {$length = $1;}
}

open FILE2, $changes or die;
my @lines = <FILE2>;
foreach my $line (@lines) {
    my @col = split (/\t/, $line);
    my $RL = length($col[3]);
    if ($col[4] !~ m/,/) {
        my $AL = length($col[4]);
        if ($col[1] <= $start) {$start = $start + ($AL-$RL);}
        if ($col[1] <= $end) {$end = $end + ($AL-$RL);}
        if ($col[1] <= $length) {$length = $length + ($AL-$RL);}
    }
    else {
        my @altbase = split (/,/,$col[4]);
        my @AD = split (/:/,$col[9]);
        my @altcoverage = split (/,/,$AD[1]);
        if ($altcoverage[2] > $altcoverage[1]) {
            my $AL = length($altbase[1]);
            if ($col[1] < $start) {$start = $start + ($AL-$RL);}
            if ($col[1] < $end) {$end = $end + ($AL-$RL);}
            if ($col[1] < $length) {$length = $length + ($AL-$RL);}
        }
        else {
            my $AL = length($altbase[0]);
            if ($col[1] < $start) {$start = $start + ($AL-$RL);}
            if ($col[1] < $end) {$end = $end + ($AL-$RL);}
            if ($col[1] < $length) {$length = $length + ($AL-$RL);}
        }
    }
}
    
print "$start\t$end\t$length\n";
EOF2
perl -e "$perlscript_IR_sizes" /path/plastid/plastid_first_check/${PI}/Solanum_plastid_assembly_2.sh /path/plastid/plastid_first_check/${PI}/${PI}.variants.vcf >${PI}/new_sizes.txt

start_of_IR=`awk '{print $1}' ${PI}/new_sizes.txt`
end_of_IR=`awk '{print $2}' ${PI}/new_sizes.txt`
start_of_IR2=`awk '{print $3}' ${PI}/new_sizes.txt`
length_of_fasta=`cat ${PI}/${PI}.plastid2.fa | tail -n -1 | awk '{print length+1}'`
###Length of LSC, IR, SSC;
LSC_length=$((start_of_IR - 1))
IR_length=$((end_of_IR - start_of_IR +1))
SSC_length=$(((start_of_IR2 - 1) - (end_of_IR +1) +1 ))
###Check if lengths were done properly
if [[ $length_of_fasta -ne $start_of_IR2 ]] ; then
    echo "quit because lengths don't match"
    exit 1
fi

###Map all reads to the final plastid sequence (by accession) using bwa
bwa index ${PI}/${ref_plastid}
bwa mem -a -M -t 8 -R ${Read_Group} ${PI}/${ref_plastid} ${fastq1} ${fastq2} >${PI}/${PI}.plastid.sam
###sort and index
samtools faidx ${PI}/${ref_plastid}
samtools view -bt ${PI}/${ref_plastid}.fai -o ${PI}/${PI}.plastid.bam ${PI}/${PI}.plastid.sam
samtools sort ${PI}/${PI}.plastid.bam ${PI}/${PI}.plastid_s
samtools view -b -q 30 -F 8 ${PI}/${PI}.plastid_s.bam >${PI}/${PI}.plastid_mpq.bam
samtools index ${PI}/${PI}.plastid_mpq.bam
rm ${PI}/${PI}.plastid.sam
rm ${PI}/${PI}.plastid.bam
rm ${PI}/${PI}.plastid_s.bam


###Use GATK to confirm that the sequence is correct
java -jar -Xmx4g /usr/local/bin/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=${PI}/${ref_plastid} OUTPUT=${PI}/${ref_plastid::(-2)}dict
java -Xmx2g -jar /usr/local/bin/picard-tools/picard.jar MarkDuplicates I=${PI}/${PI}.plastid_mpq.bam M=${PI}/metrics.${PI}.txt TMP_DIR=/dauc2/tmp MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 O=${PI}/${PI}.mpqu.bam
samtools index ${PI}/${PI}.mpqu.bam
java -jar -Xmx8g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${PI}/${ref_plastid} -I ${PI}/${PI}.mpqu.bam -o ${PI}/${PI}plastid.target_intervals.list
java -jar -Xmx8g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T IndelRealigner -R ${PI}/${ref_plastid} -I ${PI}/${PI}.mpqu.bam -targetIntervals ${PI}/${PI}plastid.target_intervals.list -o ${PI}/${PI}plastid.realigned.bam
java -jar -Xmx4g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${PI}/${ref_plastid} -I ${PI}/${PI}plastid.realigned.bam -ploidy 2 -glm BOTH -baq CALCULATE_AS_NECESSARY -dt NONE -o ${PI}/${PI}.UGvariants.vcf


###Use GATK to find coverage of plastid (LSC,1 IR, SSC)
java -jar -Xmx4g /usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar -T DepthOfCoverage -R ${PI}/${ref_plastid} -o ${PI}/${PI}_plastid_cov -I ${PI}/${PI}plastid.realigned.bam


###divide Plastid by region and find average coverage
###Then print out stats
if [ -s ${PI}plastid_final_stats ]; then #check if done; remove to redo
    rm ${PI}plastid_final_stats
fi

###separate file by region
sed -n "2,$((LSC_length+1))p" ${PI}/${PI}_plastid_cov >${PI}/${PI}_plastid_covLSC
sed -n "$((LSC_length+2)),$((LSC_length+2+IR_length))p" ${PI}/${PI}_plastid_cov >${PI}/${PI}_plastid_covIR
sed -n "$((LSC_length+3+IR_length)),$((LSC_length+3+IR_length+SSC_length))p" ${PI}/${PI}_plastid_cov >${PI}/${PI}_plastid_covSSC
	
###Find the average coverage of region
echo "Average coverage of LSC is " >>${PI}/${PI}plastid_final_stats
awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}/${PI}_plastid_covLSC >>${PI}/${PI}plastid_final_stats
echo "Average coverage of IR is " >>${PI}/${PI}plastid_final_stats
awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}/${PI}_plastid_covIR >>${PI}/${PI}plastid_final_stats
echo "Average coverage of SSC is " >>${PI}/${PI}plastid_final_stats
awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}/${PI}_plastid_covSSC >>${PI}/${PI}plastid_final_stats

###Print out regions that do not meet the Low (less than 75% (LSC, SSC or IR)) and High (1.25X (LSC, SSC, or IR)) regions
Plastid_max_SC=`awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}/${PI}_plastid_covLSC`
Pint_SC=`printf "%.0f\n" $Plastid_max_SC`
Plastid_max_IR=`awk '{sum +=$2; n++} END {if (n>0) print sum /n;}' ${PI}/${PI}_plastid_covIR`
Pint_IR=`printf "%.0f\n" $Plastid_max_IR`
MmaxSC=$((Pint_SC * 3 / 5))
PmaxSC=$((Pint_SC * 5 / 3))
MmaxIR=$((Pint_IR * 3 / 5))
PmaxIR=$((Pint_IR * 5 / 3))
cat ${PI}/${PI}_plastid_covLSC | tail -n+75 | head -n-75 | awk '$2<'$MmaxSC >>${PI}/${PI}_plastid_LowRegions
awk '$2>'$PmaxSC ${PI}/${PI}_plastid_covLSC >>${PI}/${PI}_plastid_HighRegions
cat ${PI}/${PI}_plastid_covSSC | tail -n+75 | head -n-75 | awk '$2<'$MmaxSC >>${PI}/${PI}_plastid_LowRegions
awk '$2>'$PmaxSC ${PI}/${PI}_plastid_covSSC >>${PI}/${PI}_plastid_HighRegions
cat ${PI}/${PI}_plastid_covIR | tail -n+75 | head -n-75 | awk '$2<'$MmaxIR >>${PI}/${PI}_plastid_LowRegions
awk '$2>'$PmaxIR ${PI}/${PI}_plastid_covIR >>${PI}/${PI}_plastid_HighRegions


###Find those variants that the alternate is greater than the reference
read -r -d '' perlscript_filter_plastid_vcf <<'EOF3'
use strict; use warnings;

my ($vcf, $PI) = @ARGV;

my (@col, @colon, @comma);

open FH, $vcf or die;

open redoSNP, ">", "$PI/$PI.variants.vcf";
open possibleSNP, ">", "$PI/$PI.possibleSNP.vcf";

while (my $lines = <FH>) {
  chomp $lines;
  next if ($lines =~ m/^#/); #Skip header lines with "#"
  my @col = split (/\t/, $lines); #Split the rest of the lines by tab
  my @colon = split (/:/, $col[9]); #Split the last tab by a :
  my @comma = split (/,/, $colon[1]); #Split the second variable by a ,
  if ($comma[1] > $comma[0]) {print redoSNP "$lines\n";}
  elsif (($comma[1]/$comma[0]) > 0.33) {print possibleSNP "$lines\n";}
}
EOF3
perl -e "$perlscript_filter_plastid_vcf" ${PI}/${PI}.UGvariants.vcf ${PI}

###Find any softclips not near the beginning, end or IR junctions of the sequence (wihtin 10 bases) that have a coverage greater than 30% of the coverage of the LSC
bb.softclip --bam=${PI}/${PI}plastid.realigned.bam --reference=${PI}/${PI}.plastid2.fa --outfile=${PI}/${PI}.softclip.txt

read -r -d '' perlscript_filter_softclips <<'EOF4'
use strict; use warnings;

my ($softclip, $start_of_IR, $end_of_IR, $start_of_IR2, $LSCcoverage) = @ARGV;

open FILE, $softclip or die;
while (my $lines = <FILE>) {
    chomp $lines;
    next if $. < 2;
    my @col = split(/\t/, $lines);
    if ((abs($col[1] -1)>10) && (abs($col[1] - $start_of_IR)>10) && (abs($col[1] - $end_of_IR)>10) && (abs($col[1] - $start_of_IR2)>10) && ($col[3] > (0.3 * $LSCcoverage))) {print $lines;}
}


EOF4
perl -e "$perlscript_filter_softclips" ${PI}/${PI}.softclip.txt $start_of_IR $end_of_IR $start_of_IR2 $Plastid_max_SC >${PI}/${PI}.softclip_review.txt

###Create Coverage plot vs GC content (101 size bins)
bb.coverage --reference=${PI}/${ref_plastid} --bam=${PI}/${PI}plastid.realigned.bam --bed=${PI}.plastid.bed --gc --image=${PI}.plastid_CovGC

mv *.png ${PI}/${PI}.CovGC.png

rm ${PI}/new_sizes.txt ${PI}/${PI}.mpqu.bam* ${PI}/${PI}.plastid_mpq.bam* ${PI}/metrics* ${PI}/*target_intervals.list ${PI}/${PI}_plastid_cov.sample* ${PI}/${PI}_plastid_covLSC ${PI}/${PI}_plastid_covIR ${PI}/${PI}_plastid_covSSC ${PI}/*.fa.* ${PI}/*.dict *.bed *.data
 
done ###for PI in $accessions ; do
###############################END PROGRAM###############################
